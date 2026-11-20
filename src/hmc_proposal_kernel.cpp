#include "hmc_proposal_kernel.h"
#include "factor_variables.h"
#include "index.h"
#include "transform.h"

namespace ilike
{

HMCProposalKernel::HMCProposalKernel()
:ProposalKernel(), steps(1)
{
}

HMCProposalKernel::~HMCProposalKernel()
{
}

HMCProposalKernel::HMCProposalKernel(const std::string &variable_in,
                                     const arma::mat &covariance_in,
                                     double step_in,
                                     size_t steps_in)
:ProposalKernel(), variable(variable_in), info(covariance_in, step_in*step_in), steps(steps_in)
{
  if (!(step_in > 0.0))
    Rcpp::stop("ilike::hmc: the step size is positive.");
  if (steps_in < 1)
    Rcpp::stop("ilike::hmc: at least one leapfrog step.");
}

HMCProposalKernel::HMCProposalKernel(const HMCProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void HMCProposalKernel::operator=(const HMCProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* HMCProposalKernel::duplicate() const
{
  return new HMCProposalKernel(*this);
}

ProposalKernel* HMCProposalKernel::proposal_kernel_duplicate() const
{
  return new HMCProposalKernel(*this);
}

void HMCProposalKernel::make_copy(const HMCProposalKernel &another)
{
  this->variable = another.variable;
  this->info = another.info;
  this->steps = another.steps;
  this->max_depth = another.max_depth;
  this->stan_variant = another.stan_variant;
  if (another.index != NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

arma::colvec HMCProposalKernel::gradient_at(const Particle &particle, const Parameters &base, const arma::colvec &x) const
{
  if (particle.factor_variables == NULL)
    return arma::zeros<arma::colvec>(x.n_elem);
  if (this->index == NULL)
    Rcpp::stop("ilike::hmc: the move's factors are not known.");
  // the particle's parameters with the block at x (new matrices: the
  // particle's are shared and must not change), on the model's scale
  Parameters block;
  block.set_stacked(this->variable, x, base);
  Parameters point = base;
  point.deep_overwrite_with_variables_in_argument(block);
  Parameters at = (this->transform != NULL) ? this->transform->inverse_transform(point) : point;
  std::vector<std::string> names = split_stacked_name(this->variable);
  std::vector<arma::mat> gradients = particle.factor_variables->direct_get_gradients_of_log_at(this->index.get(), names, at);
  return this->stacked_gradient(gradients, point, x.n_elem);
}

arma::colvec HMCProposalKernel::stacked_gradient(const std::vector<arma::mat> &gradients, const Parameters &point, arma::uword n) const
{
  // the model-scale gradients, stacked (their total size may differ from n
  // on the transformed scale: simplex, cov, corr, C232)
  arma::uword total = 0;
  for (const arma::mat &gradient : gradients)
    total += gradient.n_elem;
  if ((this->transform == NULL) && (total != n))
    Rcpp::stop("ilike::hmc: the gradient does not match the size of " + stacked_display_name(this->variable) + ".");
  arma::colvec stacked(total);
  arma::uword offset = 0;
  for (const arma::mat &gradient : gradients)
  {
    if (gradient.n_elem > 0)
      stacked.subvec(offset, offset + gradient.n_elem - 1) = arma::vectorise(gradient);
    offset += gradient.n_elem;
  }
  // on the transformed scale (C209): the chain rule, with log|dx/dz|
  if (this->transform != NULL)
  {
    arma::mat transformed = stacked;
    if (!this->transform->gradient_to_transformed(this->variable, point, transformed))
      Rcpp::stop("ilike::hmc on a scale given by a transform written by the user is not available "
                 "(the built-in transforms are: transform = auto, or a list of log, logit, lower, upper, simplex, cov, corr).");
    if (transformed.n_elem != n)
      Rcpp::stop("ilike::hmc: the gradient does not match the size of " + stacked_display_name(this->variable) + ".");
    return arma::vectorise(transformed);
  }
  return stacked;
}

Parameters HMCProposalKernel::trajectory(RandomNumberGenerator &rng,
                                         const Particle &particle,
                                         double &log_ratio) const
{
  // (on the kernel's scale: transformed, with a transform, C209)
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x0 = current.get_colvec(this->variable);
  arma::uword n = x0.n_elem;
  this->check_size(n);
  double eps = this->step_size();
  size_t number_of_steps = this->steps;
  // a particle's own (step, steps) (esjd in SMC, C211)
  if (const double* tuning = this->tuning_values(particle))
  {
    eps = tuning[0];
    number_of_steps = size_t(std::max(1.0, std::round(tuning[1])));
  }
  arma::colvec r = this->draw_momentum(rng, n);
  auto apply_c = [&](const arma::colvec &v) -> arma::colvec { return this->velocity(v); };
  auto kinetic = [&](const arma::colvec &v) -> double { return this->kinetic(v); };
  double start = kinetic(r);

  arma::colvec x = x0;
  arma::colvec g = this->gradient_at(particle, current, x);
  bool finite = g.is_finite();
  for (size_t l = 0; (l < number_of_steps) && finite; ++l)
  {
    r = r + (0.5*eps)*g;
    x = x + eps*apply_c(r);
    g = this->gradient_at(particle, current, x);
    r = r + (0.5*eps)*g;
    finite = x.is_finite() && g.is_finite() && r.is_finite();
  }

  Parameters output;
  if (!finite)
  {
    // (the trajectory reached a point of density 0, or diverged: no move)
    log_ratio = -arma::datum::inf;
    output.set_stacked(this->variable, x0, current);
    return output;
  }
  log_ratio = start - kinetic(r);
  output.set_stacked(this->variable, x, current);
  return output;
}

void HMCProposalKernel::check_size(arma::uword n) const
{
  const arma::mat &covariance = this->info.get_covariance();
  if (covariance.n_rows != 1)
    check_proposal_size(this->variable, n, covariance.n_rows, "mass (covariance)");
}

// r ~ N(0, C^-1): r = L^-T z with C = L L' (c^-1/2 z for C = c I)
arma::colvec HMCProposalKernel::draw_momentum(RandomNumberGenerator &rng, arma::uword n) const
{
  const arma::mat &covariance = this->info.get_covariance();
  arma::colvec z = arma::vectorise(rnorm(rng, std::pair<size_t, size_t>(n, 1)));
  if (covariance.n_rows == 1)
    return z/std::sqrt(covariance(0, 0));
  return arma::solve(arma::trimatu(this->info.get_chol_lower().t()), z);
}

arma::colvec HMCProposalKernel::velocity(const arma::colvec &r) const
{
  const arma::mat &covariance = this->info.get_covariance();
  return (covariance.n_rows == 1) ? arma::colvec(covariance(0, 0)*r) : arma::colvec(covariance*r);
}

double HMCProposalKernel::kinetic(const arma::colvec &r) const
{
  return 0.5*arma::dot(r, this->velocity(r));
}

double HMCProposalKernel::log_density_at(const Particle &particle,
                                         const Parameters &base,
                                         const arma::colvec &x,
                                         const Index* index,
                                         const Index* new_outputs_index,
                                         std::shared_ptr<Particle> &point) const
{
  point.reset();
  if (!x.is_finite())
    return -arma::datum::inf;
  Parameters block;
  block.set_stacked(this->variable, x, base);
  point = std::make_shared<Particle>(this->proposed_from(particle, block, new_outputs_index));
  double log_density = point->evaluate_likelihoods(index);
  if (this->transform != NULL)
  {
    Parameters at = base;
    at.deep_overwrite_with_variables_in_argument(block);
    log_density = log_density + this->transform->log_abs_inverse_jacobian_determinant(at);
  }
  return log_density;
}

double HMCProposalKernel::log_density_and_gradient_at(const Particle &particle,
                                                       const Parameters &base,
                                                       const arma::colvec &x,
                                                       const Index* index,
                                                       const Index* new_outputs_index,
                                                       std::shared_ptr<Particle> &point,
                                                       arma::colvec &gradient) const
{
  point.reset();
  if (!x.is_finite() || (particle.factor_variables == NULL))
  {
    gradient = arma::zeros<arma::colvec>(x.n_elem);
    return x.is_finite() ? this->log_density_at(particle, base, x, index, new_outputs_index, point) : -arma::datum::inf;
  }
  Parameters block;
  block.set_stacked(this->variable, x, base);
  point = std::make_shared<Particle>(this->proposed_from(particle, block, new_outputs_index));
  std::vector<std::string> names = split_stacked_name(this->variable);
  std::vector<arma::mat> gradients;
  double log_density = point->factor_variables->evaluate_likelihoods_with_gradients(index, names, gradients);
  point->target_evaluated = log_density;
  Parameters at = base;
  at.deep_overwrite_with_variables_in_argument(block);
  if (this->transform != NULL)
    log_density = log_density + this->transform->log_abs_inverse_jacobian_determinant(at);
  if (std::isfinite(log_density))
    gradient = this->stacked_gradient(gradients, at, x.n_elem);
  else
    gradient = arma::zeros<arma::colvec>(x.n_elem);
  return log_density;
}

double HMCProposalKernel::target_at(const Particle &particle, const Parameters &base, const arma::colvec &x) const
{
  if (this->index == NULL)
    Rcpp::stop("ilike::hmc: the move's factors are not known.");
  std::shared_ptr<Particle> point;
  return this->log_density_at(particle, base, x, this->index.get(), this->index.get(), point);
}

double HMCProposalKernel::log_jacobian(const Parameters &base) const
{
  return (this->transform != NULL) ? this->transform->log_abs_inverse_jacobian_determinant(base) : 0.0;
}

Particle HMCProposalKernel::move(RandomNumberGenerator &rng,
                                 const Particle &particle,
                                 const Index* new_outputs_index) const
{
  double log_ratio;
  Parameters proposed = this->trajectory(rng, particle, log_ratio);
  Particle proposed_particle = this->proposed_from(particle, proposed, new_outputs_index);
  proposed_particle.proposal_serial = particle.proposal_serial + 1;
  proposed_particle.proposal_log_ratio = log_ratio;
  return proposed_particle;
}

void HMCProposalKernel::move_into(RandomNumberGenerator &rng,
                                  const Particle &particle,
                                  const Index* new_outputs_index,
                                  Particle &scratch) const
{
  double log_ratio;
  Parameters proposed = this->trajectory(rng, particle, log_ratio);
  this->into_from(particle, proposed, new_outputs_index, scratch);
  scratch.proposal_serial = particle.proposal_serial + 1;
  scratch.proposal_log_ratio = log_ratio;
}

Parameters HMCProposalKernel::simulate(RandomNumberGenerator &rng,
                                       const Particle &particle) const
{
  double log_ratio;
  return this->trajectory(rng, particle, log_ratio);
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const Particle &particle) const
{
  Rcpp::stop("ilike::hmc with subsampling is not available.");
}

Parameters HMCProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                 const std::string &variable,
                                                 const Particle &particle) const
{
  Rcpp::stop("ilike::hmc with subsampling is not available.");
}

// log q(proposed | old), for the Metropolis-Hastings ratio: of the two
// particles, the one made by this proposal from the other has a serial one
// more; log q(new | old) is taken as 0 and log q(old | new) as K(r) - K(r*),
// so that the ratio's proposal term is exp(K(r) - K(r*)).
double HMCProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                   const Particle &old_particle) const
{
  if (proposed_particle.proposal_serial == old_particle.proposal_serial + 1)
    return 0.0;
  if (old_particle.proposal_serial == proposed_particle.proposal_serial + 1)
    return old_particle.proposal_log_ratio;
  Rcpp::stop("ilike::hmc: the density of the HMC proposal between two arbitrary points is not available "
             "(hmc is used in metropolis_hastings and unadjusted moves).");
}

double HMCProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                             const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

arma::mat HMCProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                      const Particle &proposed_particle,
                                                      const Particle &old_particle)
{
  Rcpp::stop("ilike::hmc: the gradient of the proposal's density is not available.");
}

arma::mat HMCProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                const Particle &proposed_particle,
                                                                const Particle &old_particle)
{
  Rcpp::stop("ilike::hmc: the gradient of the proposal's density is not available.");
}

void HMCProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* HMCProposalKernel::simulate_gradient_estimator_output() const
{
  // (the gradients are taken along the trajectory, not kept with the particle)
  return NULL;
}

std::vector<const ProposalKernel*> HMCProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

// a private copy of the index (the caller owns and deletes the one given)
void HMCProposalKernel::set_index(Index* index_in)
{
  if (index_in != NULL)
    this->index.reset(index_in->duplicate());
  else
    this->index.reset();
}

void HMCProposalKernel::set_index_if_null(Index* index_in)
{
  if ((this->index == NULL) && (index_in != NULL))
    this->index.reset(index_in->duplicate());
}

bool HMCProposalKernel::can_be_evaluated() const
{
  return true;
}

void HMCProposalKernel::set_data(Data* data_in)
{
}

GaussianProposalInfo* HMCProposalKernel::get_gaussian_proposal_info(const std::string &variable_in)
{
  return (variable_in == this->variable) ? &this->info : NULL;
}

}
