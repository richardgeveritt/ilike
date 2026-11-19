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
    next.log_density = kernel->log_density_at(particle, base, next.x, this->index.get(),
                                              this->get_new_outputs_index(this->index.get()), next.particle);
    if (std::isfinite(next.log_density))
      next.gradient = kernel->gradient_at(particle, base, next.x);
    else
      next.gradient = arma::zeros<arma::colvec>(next.x.n_elem);
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
    tree.n = total;
  }
  return tree;
}

Particle NoUTurnMCMC::transition(RandomNumberGenerator &rng,
                                 const Particle &particle,
                                 double &statistic,
                                 size_t &depth,
                                 bool &divergent) const
{
  const HMCProposalKernel* kernel = this->hamiltonian();
  const Parameters &base = particle.get_transformed_parameters(kernel);
  State start;
  start.x = base.get_colvec(kernel->block());
  kernel->check_size(start.x.n_elem);
  double eps = kernel->step_size();
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
  return this->transition(rng, particle, statistic, depth, divergent);
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
  this->proposal->prepare_mcmc_adaptation();
  double statistic;
  size_t depth;
  bool divergent;
  Particle next = this->transition(rng, particle, statistic, depth, divergent);
  scratch = particle;
  particle = next;
  AdaptationStep step;
  step.previous = &scratch;
  step.proposed = &particle;
  step.current = &particle;
  step.acceptance_probability = statistic;
  this->proposal->mcmc_adapt(step);
}

Particle NoUTurnMCMC::subsample_move(RandomNumberGenerator &rng,
                                     const Particle &particle) const
{
  Rcpp::stop("nuts with subsampling is not available.");
}

}
