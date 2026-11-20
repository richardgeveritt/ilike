#include "no_u_turn_mcmc.h"
#include "hmc_proposal_kernel.h"
#include "mcmc_adaptor.h"
#include "index.h"

namespace ilike
{

NoUTurnMCMC::NoUTurnMCMC()
:MetropolisHastingsMCMC()
{
}

NoUTurnMCMC::NoUTurnMCMC(size_t number_of_iterations_in,
                         HMCProposalKernel* proposal_in)
:MetropolisHastingsMCMC(number_of_iterations_in, proposal_in)
{
}

NoUTurnMCMC::NoUTurnMCMC(MCMCTermination* termination_in,
                         HMCProposalKernel* proposal_in)
:MetropolisHastingsMCMC(termination_in, proposal_in)
{
}

NoUTurnMCMC::~NoUTurnMCMC()
{
}

NoUTurnMCMC::NoUTurnMCMC(const NoUTurnMCMC &another)
:MetropolisHastingsMCMC(another)
{
}

void NoUTurnMCMC::operator=(const NoUTurnMCMC &another)
{
  if (this == &another)
    return;
  MetropolisHastingsMCMC::operator=(another);
}

Kernel* NoUTurnMCMC::duplicate() const
{
  return new NoUTurnMCMC(*this);
}

MCMC* NoUTurnMCMC::mcmc_duplicate() const
{
  return new NoUTurnMCMC(*this);
}

MetropolisHastingsMCMC* NoUTurnMCMC::metropolis_hastings_mcmc_duplicate() const
{
  return new NoUTurnMCMC(*this);
}

const HMCProposalKernel* NoUTurnMCMC::hamiltonian() const
{
  const HMCProposalKernel* kernel = dynamic_cast<const HMCProposalKernel*>(this->proposal.get());
  if (kernel == NULL)
    Rcpp::stop("nuts: the move has no Hamiltonian system.");
  return kernel;
}

bool NoUTurnMCMC::no_u_turn(const State &minus, const State &plus) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  arma::colvec span = plus.x - minus.x;
  return (arma::dot(span, kernel->velocity(minus.r)) >= 0.0) &&
         (arma::dot(span, kernel->velocity(plus.r)) >= 0.0);
}

NoUTurnMCMC::Tree NoUTurnMCMC::build_tree(RandomNumberGenerator &rng,
                                          const Particle &particle,
                                          const Parameters &base,
                                          const State &from,
                                          double log_u,
                                          int direction,
                                          size_t depth,
                                          double eps,
                                          double joint0) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  if (depth == 0)
  {
    // one leapfrog step in the direction given
    double step = double(direction)*eps;
    State next;
    next.r = from.r + (0.5*step)*from.gradient;
    next.x = from.x + step*kernel->velocity(next.r);
    // the log density and its gradient in one pass over the terms (C222)
    next.log_density = kernel->log_density_and_gradient_at(particle, base, next.x, this->index.get(),
                                                           this->get_new_outputs_index(this->index.get()), next.particle,
                                                           next.gradient);
    next.r = next.r + (0.5*step)*next.gradient;
    double joint = next.log_density - kernel->kinetic(next.r);
    if (!next.gradient.is_finite() || !next.r.is_finite())
      joint = -arma::datum::inf;
    Tree tree;
    tree.minus = next;
    tree.plus = next;
    tree.candidate = next;
    tree.n = (log_u <= joint) ? 1.0 : 0.0;
    tree.s = (log_u < 1000.0 + joint);
    tree.divergent = !tree.s;
    double change = joint - joint0;
    tree.alpha = std::isnan(change) ? 0.0 : std::min(1.0, std::exp(change));
    tree.n_alpha = 1.0;
    tree.steps = 1;
    return tree;
  }
  // the left and right subtrees, built implicitly
  Tree tree = this->build_tree(rng, particle, base, from, log_u, direction, depth - 1, eps, joint0);
  if (tree.s)
  {
    Tree second = this->build_tree(rng, particle, base, (direction == -1) ? tree.minus : tree.plus,
                                   log_u, direction, depth - 1, eps, joint0);
    if (direction == -1)
      tree.minus = second.minus;
    else
      tree.plus = second.plus;
    double total = tree.n + second.n;
    if ((total > 0.0) && (runif(rng) < second.n/total))
      tree.candidate = second.candidate;
    tree.alpha = tree.alpha + second.alpha;
    tree.n_alpha = tree.n_alpha + second.n_alpha;
    tree.s = second.s && this->no_u_turn(tree.minus, tree.plus);
    tree.divergent = tree.divergent || second.divergent;
    tree.steps = tree.steps + second.steps;
    tree.n = total;
  }
  return tree;
}

// --- Stan's variant (C221) ----------------------------------------------------

static double log_sum_exp(double a, double b)
{
  if (a == -arma::datum::inf)
    return b;
  if (b == -arma::datum::inf)
    return a;
  double m = std::max(a, b);
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

bool NoUTurnMCMC::criterion(const arma::colvec &r_minus, const arma::colvec &r_plus, const arma::colvec &rho) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  return (arma::dot(kernel->velocity(r_minus), rho) > 0.0) && (arma::dot(kernel->velocity(r_plus), rho) > 0.0);
}

bool NoUTurnMCMC::merge_criterion(const State &left_minus, const State &left_plus, const arma::colvec &left_rho,
                                  const State &right_minus, const State &right_plus, const arma::colvec &right_rho) const
{
  // the whole trajectory, and (Stan, since 2.26) the left trajectory with
  // the first state of the right one, and the right with the last of the left
  return this->criterion(left_minus.r, right_plus.r, left_rho + right_rho) &&
         this->criterion(left_minus.r, right_minus.r, left_rho + right_minus.r) &&
         this->criterion(left_plus.r, right_plus.r, right_rho + left_plus.r);
}

NoUTurnMCMC::WeightedTree NoUTurnMCMC::build_weighted_tree(RandomNumberGenerator &rng,
                                                           const Particle &particle,
                                                           const Parameters &base,
                                                           const State &from,
                                                           int direction,
                                                           size_t depth,
                                                           double eps,
                                                           double joint0) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  if (depth == 0)
  {
    // one leapfrog step; the state's weight exp(H0 - H) = exp(joint - joint0)
    double step = double(direction)*eps;
    State next;
    next.r = from.r + (0.5*step)*from.gradient;
    next.x = from.x + step*kernel->velocity(next.r);
    // the log density and its gradient in one pass over the terms (C222)
    next.log_density = kernel->log_density_and_gradient_at(particle, base, next.x, this->index.get(),
                                                           this->get_new_outputs_index(this->index.get()), next.particle,
                                                           next.gradient);
    next.r = next.r + (0.5*step)*next.gradient;
    double joint = next.log_density - kernel->kinetic(next.r);
    if (!next.gradient.is_finite() || !next.r.is_finite() || std::isnan(joint))
      joint = -arma::datum::inf;
    WeightedTree tree;
    tree.minus = next;
    tree.plus = next;
    tree.candidate = next;
    tree.log_w = joint - joint0;
    tree.rho = next.r;
    tree.divergent = !(joint0 - joint <= 1000.0); // (H - H0 above 1000, or not finite)
    tree.s = !tree.divergent;
    tree.sum_alpha = std::isfinite(joint) ? std::min(1.0, std::exp(joint - joint0)) : 0.0;
    tree.steps = 1;
    return tree;
  }
  WeightedTree first = this->build_weighted_tree(rng, particle, base, from, direction, depth - 1, eps, joint0);
  if (!first.s)
    return first;
  WeightedTree second = this->build_weighted_tree(rng, particle, base, (direction == -1) ? first.minus : first.plus,
                                                  direction, depth - 1, eps, joint0);
  WeightedTree tree = first;
  tree.steps = first.steps + second.steps;
  tree.sum_alpha = first.sum_alpha + second.sum_alpha;
  tree.divergent = first.divergent || second.divergent;
  if (!second.s)
  {
    tree.s = false;
    return tree;
  }
  if (direction == -1)
    tree.minus = second.minus;
  else
    tree.plus = second.plus;
  // the candidate in proportion to the weights (uniform progressive sampling)
  double log_w = log_sum_exp(first.log_w, second.log_w);
  if ((second.log_w > -arma::datum::inf) && (std::log(runif(rng)) < second.log_w - log_w))
    tree.candidate = second.candidate;
  tree.log_w = log_w;
  tree.rho = first.rho + second.rho;
  // left: the earlier of the two in time
  const WeightedTree &left = (direction == -1) ? second : first;
  const WeightedTree &right = (direction == -1) ? first : second;
  tree.s = this->merge_criterion(left.minus, left.plus, left.rho, right.minus, right.plus, right.rho);
  return tree;
}

Particle NoUTurnMCMC::stan_transition(RandomNumberGenerator &rng,
                                      const Particle &particle,
                                      double &statistic,
                                      size_t &depth,
                                      bool &divergent,
                                      size_t &leapfrog_steps) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  const Parameters &base = particle.get_transformed_parameters(kernel);
  State start;
  start.x = base.get_colvec(kernel->block());
  kernel->check_size(start.x.n_elem);
  double eps = kernel->step_size();
  if (const double* tuning = kernel->tuning_values(particle))
    eps = tuning[0];
  start.r = kernel->draw_momentum(rng, start.x.n_elem);
  start.log_density = particle.get_cached_likelihoods(this->index.get()) + kernel->log_jacobian(base);
  start.gradient = kernel->gradient_at(particle, base, start.x);
  double joint0 = start.log_density - kernel->kinetic(start.r);

  State minus = start, plus = start, chosen = start;
  arma::colvec rho = start.r;
  double log_w = 0.0; // (the start's weight, exp(0))
  double sum_alpha = 0.0;
  bool s = std::isfinite(joint0) && start.gradient.is_finite();
  statistic = 0.0;
  depth = 0;
  divergent = false;
  leapfrog_steps = 0;
  size_t max_depth = kernel->get_max_depth();
  while (s && (depth < max_depth))
  {
    int direction = (runif(rng) < 0.5) ? -1 : 1;
    WeightedTree tree = this->build_weighted_tree(rng, particle, base, (direction == -1) ? minus : plus,
                                                  direction, depth, eps, joint0);
    leapfrog_steps = leapfrog_steps + tree.steps;
    sum_alpha = sum_alpha + tree.sum_alpha;
    divergent = divergent || tree.divergent;
    depth = depth + 1;
    if (!tree.s)
      break;
    // biased progressive sampling: the new subtree's candidate with
    // probability min(1, w_new/w_old), which favours states far from the start
    if ((tree.log_w > -arma::datum::inf) && (std::log(runif(rng)) < tree.log_w - log_w))
      chosen = tree.candidate;
    log_w = log_sum_exp(log_w, tree.log_w);
    const State old_minus = minus, old_plus = plus;
    const arma::colvec old_rho = rho;
    if (direction == -1)
    {
      minus = tree.minus;
      s = this->merge_criterion(tree.minus, tree.plus, tree.rho, old_minus, old_plus, old_rho);
    }
    else
    {
      plus = tree.plus;
      s = this->merge_criterion(old_minus, old_plus, old_rho, tree.minus, tree.plus, tree.rho);
    }
    rho = old_rho + tree.rho;
  }
  // the acceptance statistic: the average over every state built (as Stan)
  statistic = (leapfrog_steps > 0) ? sum_alpha/double(leapfrog_steps) : 0.0;

  if (chosen.particle == NULL)
    return particle; // (the current point)
  Particle output = *chosen.particle;
  output.set_acceptance(this->proposal.get(), true);
  this->refresh_cut_factors(output);
  this->check_carried_over_factors(output, this->index.get());
  return output;
}

Particle NoUTurnMCMC::transition(RandomNumberGenerator &rng,
                                 const Particle &particle,
                                 double &statistic,
                                 size_t &depth,
                                 bool &divergent,
                                 size_t &leapfrog_steps) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  if (kernel->get_stan_variant())
    return this->stan_transition(rng, particle, statistic, depth, divergent, leapfrog_steps);
  const Parameters &base = particle.get_transformed_parameters(kernel);
  State start;
  start.x = base.get_colvec(kernel->block());
  kernel->check_size(start.x.n_elem);
  double eps = kernel->step_size();
  // a particle's own step (esjd in SMC, C211)
  if (const double* tuning = kernel->tuning_values(particle))
    eps = tuning[0];
  start.r = kernel->draw_momentum(rng, start.x.n_elem);
  start.log_density = particle.get_cached_likelihoods(this->index.get()) + kernel->log_jacobian(base);
  start.gradient = kernel->gradient_at(particle, base, start.x);
  double joint0 = start.log_density - kernel->kinetic(start.r);
  // the slice: u ~ U(0, exp(joint0))
  double log_u = joint0 + std::log(runif(rng));

  State minus = start, plus = start, chosen = start;
  double n = 1.0;
  bool s = std::isfinite(joint0) && start.gradient.is_finite();
  statistic = 0.0;
  depth = 0;
  divergent = false;
  leapfrog_steps = 0;
  size_t max_depth = kernel->get_max_depth();
  while (s && (depth < max_depth))
  {
    int direction = (runif(rng) < 0.5) ? -1 : 1;
    Tree tree = (direction == -1)
      ? this->build_tree(rng, particle, base, minus, log_u, direction, depth, eps, joint0)
      : this->build_tree(rng, particle, base, plus, log_u, direction, depth, eps, joint0);
    if (direction == -1)
      minus = tree.minus;
    else
      plus = tree.plus;
    if (tree.s && (runif(rng) < std::min(1.0, tree.n/n)))
      chosen = tree.candidate;
    n = n + tree.n;
    s = tree.s && this->no_u_turn(minus, plus);
    statistic = (tree.n_alpha > 0.0) ? tree.alpha/tree.n_alpha : 0.0;
    divergent = divergent || tree.divergent;
    leapfrog_steps = leapfrog_steps + tree.steps;
    depth = depth + 1;
  }

  if (chosen.particle == NULL)
    return particle; // (the current point)
  Particle output = *chosen.particle;
  output.set_acceptance(this->proposal.get(), true);
  this->refresh_cut_factors(output);
  this->check_carried_over_factors(output, this->index.get());
  return output;
}

Particle NoUTurnMCMC::move(RandomNumberGenerator &rng,
                           const Particle &particle) const
{
  double statistic;
  size_t depth;
  bool divergent;
  size_t steps;
  return this->transition(rng, particle, statistic, depth, divergent, steps);
}

void NoUTurnMCMC::move_in_place(RandomNumberGenerator &rng,
                                Particle &particle) const
{
  particle = this->move(rng, particle);
}

void NoUTurnMCMC::move_in_place(RandomNumberGenerator &rng,
                                Particle &particle,
                                Particle &scratch) const
{
  particle = this->move(rng, particle);
}

void NoUTurnMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                         Particle &particle,
                                         Particle &scratch)
{
  // the adaptation step with Algorithm 6's statistic as the acceptance
  // probability (the adaptors read alpha, and the states before and after)
  double row[4];
  this->diagnostic_move_in_place(rng, particle, scratch, true, row);
}

void NoUTurnMCMC::diagnostic_names(std::vector<std::string> &names) const
{
  names = {"tree_depth", "leapfrog_steps", "divergent", "accept_stat"};
}

bool NoUTurnMCMC::warming_up() const
{
  return this->proposal->mcmc_warming_up();
}

void NoUTurnMCMC::diagnostic_move_in_place(RandomNumberGenerator &rng,
                                           Particle &particle,
                                           Particle &scratch,
                                           bool adapting,
                                           double* row)
{
  // (the same random numbers as move_in_place / adaptive_move_in_place)
  if (adapting)
    this->proposal->prepare_mcmc_adaptation();
  double statistic;
  size_t depth;
  bool divergent;
  size_t steps;
  Particle next = this->transition(rng, particle, statistic, depth, divergent, steps);
  row[0] = double(depth);
  row[1] = double(steps);
  row[2] = divergent ? 1.0 : 0.0;
  row[3] = statistic;
  if (!adapting)
  {
    particle = next;
    return;
  }
  scratch = particle;
  particle = next;
  AdaptationStep step;
  step.previous = &scratch;
  step.proposed = &particle;
  step.current = &particle;
  step.acceptance_probability = statistic;
  step.rng = &rng;
  this->proposal->mcmc_adapt(step);
}

Particle NoUTurnMCMC::subsample_move(RandomNumberGenerator &rng,
                                     const Particle &particle) const
{
  Rcpp::stop("nuts with subsampling is not available.");
}

}
