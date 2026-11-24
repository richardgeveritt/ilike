#include "pcn_langevin_proposal_kernel.h"
#include "gradient_estimator.h"
#include "gradient_estimator_output.h"
#include "particle.h"
#include "distributions.h"

namespace ilike
{

PCNLangevinProposalKernel::PCNLangevinProposalKernel()
:LangevinProposalKernel()
{
}

PCNLangevinProposalKernel::~PCNLangevinProposalKernel()
{
}

PCNLangevinProposalKernel::PCNLangevinProposalKernel(const std::string &variable_in, const arma::colvec &mean_in,
                                                     const arma::mat &covariance_in, double delta_in,
                                                     GradientEstimator* gradient_estimator_in)
:LangevinProposalKernel(variable_in, covariance_in, gradient_estimator_in), variable(variable_in), mean(mean_in),
 covariance(covariance_in)
{
  if (!(delta_in > 0.0))
    Rcpp::stop("ilike::pcn_mala: delta must be positive.");
  if ((covariance_in.n_rows != mean_in.n_elem) || (covariance_in.n_cols != mean_in.n_elem))
    Rcpp::stop("ilike::pcn_mala: the covariance must be a square matrix with as many rows as the mean has elements.");
  if (!arma::chol(this->lower_chol, covariance_in, "lower"))
    Rcpp::stop("ilike::pcn_mala: the covariance is not positive definite.");
  this->inverse_covariance = arma::inv_sympd(covariance_in);
  this->a = (2.0 - delta_in)/(2.0 + delta_in);
  this->b = 2.0*delta_in/(2.0 + delta_in);
  this->c = std::sqrt(8.0*delta_in)/(2.0 + delta_in);
  this->inverse_step_covariance = this->inverse_covariance/(this->c*this->c);
  this->log_det_step_covariance = arma::log_det_sympd(covariance_in) + double(mean_in.n_elem)*std::log(this->c*this->c);
}

PCNLangevinProposalKernel::PCNLangevinProposalKernel(const PCNLangevinProposalKernel &another)
:LangevinProposalKernel(another)
{
  this->make_copy(another);
}

void PCNLangevinProposalKernel::operator=(const PCNLangevinProposalKernel &another)
{
  if (this == &another)
    return;
  LangevinProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* PCNLangevinProposalKernel::duplicate() const
{
  return new PCNLangevinProposalKernel(*this);
}

ProposalKernel* PCNLangevinProposalKernel::proposal_kernel_duplicate() const
{
  return new PCNLangevinProposalKernel(*this);
}

void PCNLangevinProposalKernel::make_copy(const PCNLangevinProposalKernel &another)
{
  this->variable = another.variable;
  this->mean = another.mean;
  this->covariance = another.covariance;
  this->inverse_covariance = another.inverse_covariance;
  this->lower_chol = another.lower_chol;
  this->a = another.a;
  this->b = another.b;
  this->c = another.c;
  this->inverse_step_covariance = another.inverse_step_covariance;
  this->log_det_step_covariance = another.log_det_step_covariance;
}

arma::colvec PCNLangevinProposalKernel::proposal_mean(const Particle &particle) const
{
  arma::colvec u = particle.get_transformed_parameters(this).get_colvec(this->variable);
  if (u.n_elem != this->mean.n_elem)
    Rcpp::stop("ilike::pcn_mala: the variable has " + std::to_string(u.n_elem) + " element(s), the mean " +
               std::to_string(this->mean.n_elem) + ".");
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  arma::colvec gradient = arma::vectorise(gradient_estimator_output->get_gradient_of_log(this->variable, this->index.get(), particle));
  // DPhi = -(grad log pi - grad log N(u; m, C)) = -grad log pi - C^-1 (u - m)
  arma::colvec dphi = -gradient - this->inverse_covariance*(u - this->mean);
  return this->mean + this->a*(u - this->mean) - this->b*(this->covariance*dphi);
}

double PCNLangevinProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                           const Particle &old_particle) const
{
  return dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(this->variable),
                               this->proposal_mean(old_particle), this->inverse_step_covariance,
                               this->log_det_step_covariance);
}

Parameters PCNLangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  Parameters output;
  output.set_stacked(this->variable, rmvnorm_using_lower_chol(rng, this->proposal_mean(particle), this->lower_chol, this->c),
                     particle.get_transformed_parameters(this));
  return output;
}

void PCNLangevinProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                              const Particle &particle,
                                              ProposalTarget &target) const
{
  target.write(this->variable, rmvnorm_using_lower_chol(rng, this->proposal_mean(particle), this->lower_chol, this->c),
               particle.get_transformed_parameters(this));
}

}
