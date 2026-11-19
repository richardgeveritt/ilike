#include "hmc_proposal_kernel.h"
#include "factor_variables.h"
#include "index.h"

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
  if (another.index != NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

arma::colvec HMCProposalKernel::gradient_at(const Particle &particle, const arma::colvec &x) const
{
  if (particle.factor_variables == NULL)
    return arma::zeros<arma::colvec>(x.n_elem);
  if (this->index == NULL)
    Rcpp::stop("ilike::hmc: the move's factors are not known.");
  // the particle's parameters with the block at x (new matrices: the
  // particle's are shared and must not change)
  Parameters block;
  block.set_stacked(this->variable, x, particle.parameters);
  Parameters at = particle.parameters;
  at.deep_overwrite_with_variables_in_argument(block);
  std::vector<std::string> names = split_stacked_name(this->variable);
  std::vector<arma::mat> gradients = particle.factor_variables->direct_get_gradients_of_log_at(this->index.get(), names, at);
  arma::colvec stacked(x.n_elem);
  arma::uword offset = 0;
  for (const arma::mat &gradient : gradients)
  {
    if (offset + gradient.n_elem > x.n_elem)
      Rcpp::stop("ilike::hmc: the gradient does not match the size of " + stacked_display_name(this->variable) + ".");
    if (gradient.n_elem > 0)
      stacked.subvec(offset, offset + gradient.n_elem - 1) = arma::vectorise(gradient);
    offset += gradient.n_elem;
  }
  return stacked;
}

Parameters HMCProposalKernel::trajectory(RandomNumberGenerator &rng,
                                         const Particle &particle,
                                         double &log_ratio) const
{
  if (this->transform != NULL)
    Rcpp::stop("ilike::hmc on a transformed space is not available yet.");
  const Parameters &current = particle.parameters;
  arma::colvec x0 = current.get_colvec(this->variable);
  arma::uword n = x0.n_elem;
  const arma::mat &covariance = this->info.get_covariance();
  bool scalar = (covariance.n_rows == 1);
  if (!scalar)
    check_proposal_size(this->variable, n, covariance.n_rows, "mass (covariance)");
  double eps = std::sqrt(this->info.get_double_scale());

  // r ~ N(0, C^-1): r = L^-T z with C = L L' (c^-1/2 z for C = c I)
  arma::colvec z = arma::vectorise(rnorm(rng, std::pair<size_t, size_t>(n, 1)));
  arma::colvec r;
  if (scalar)
    r = z/std::sqrt(covariance(0, 0));
  else
    r = arma::solve(arma::trimatu(this->info.get_chol_lower().t()), z);
  auto apply_c = [&](const arma::colvec &v) -> arma::colvec {
    return scalar ? arma::colvec(covariance(0, 0)*v) : arma::colvec(covariance*v);
  };
  auto kinetic = [&](const arma::colvec &v) -> double { return 0.5*arma::dot(v, apply_c(v)); };
  double start = kinetic(r);

  arma::colvec x = x0;
  arma::colvec g = this->gradient_at(particle, x);
  bool finite = g.is_finite();
  for (size_t l = 0; (l < this->steps) && finite; ++l)
  {
    r = r + (0.5*eps)*g;
    x = x + eps*apply_c(r);
    g = this->gradient_at(particle, x);
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
