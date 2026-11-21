#include "pcn_proposal_kernel.h"

namespace ilike
{

PCNProposalKernel::PCNProposalKernel()
:ProposalKernel(), rho(0.0), log_det_step_cov(0.0)
{
}

PCNProposalKernel::~PCNProposalKernel()
{
}

PCNProposalKernel::PCNProposalKernel(const std::string &variable_in, const arma::colvec &mean_in,
                                     const arma::mat &covariance_in, double rho_in)
:ProposalKernel(), variable(variable_in), mean(mean_in), rho(rho_in)
{
  if (!(rho_in >= 0.0) || !(rho_in < 1.0))
    Rcpp::stop("ilike::pcn: rho must be in [0, 1).");
  if ((covariance_in.n_rows != mean_in.n_elem) || (covariance_in.n_cols != mean_in.n_elem))
    Rcpp::stop("ilike::pcn: the covariance must be a square matrix with as many rows as the mean has elements.");
  if (!arma::chol(this->lower_chol, covariance_in, "lower"))
    Rcpp::stop("ilike::pcn: the covariance is not positive definite.");
  arma::mat step_cov = (1.0 - rho_in*rho_in)*covariance_in;
  this->inverse_step_cov = arma::inv_sympd(step_cov);
  this->log_det_step_cov = arma::log_det_sympd(step_cov);
}

PCNProposalKernel::PCNProposalKernel(const PCNProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void PCNProposalKernel::operator=(const PCNProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* PCNProposalKernel::duplicate() const
{
  return new PCNProposalKernel(*this);
}

ProposalKernel* PCNProposalKernel::proposal_kernel_duplicate() const
{
  return new PCNProposalKernel(*this);
}

void PCNProposalKernel::make_copy(const PCNProposalKernel &another)
{
  this->variable = another.variable;
  this->mean = another.mean;
  this->rho = another.rho;
  this->lower_chol = another.lower_chol;
  this->inverse_step_cov = another.inverse_step_cov;
  this->log_det_step_cov = another.log_det_step_cov;
}

double PCNProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                         const Particle &old_particle) const
{
  // log q(proposed | old) = log N(proposed; m + rho (old - m), (1 - rho^2) cov)
  arma::colvec x = proposed_particle.get_transformed_parameters(this).get_colvec(this->variable);
  arma::colvec u = old_particle.get_transformed_parameters(this).get_colvec(this->variable);
  if ((x.n_elem != this->mean.n_elem) || (u.n_elem != this->mean.n_elem))
    Rcpp::stop("ilike::pcn: the variable has " + std::to_string(x.n_elem) + " element(s), the mean " +
               std::to_string(this->mean.n_elem) + ".");
  return dmvnorm_using_precomp(x, this->mean + this->rho*(u - this->mean), this->inverse_step_cov, this->log_det_step_cov);
}

double PCNProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                   const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters PCNProposalKernel::simulate(RandomNumberGenerator &rng,
                                                             const Particle &particle) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec u = current.get_colvec(this->variable);
  if (u.n_elem != this->mean.n_elem)
    Rcpp::stop("ilike::pcn: the variable has " + std::to_string(u.n_elem) + " element(s), the mean " +
               std::to_string(this->mean.n_elem) + ".");
  arma::colvec x = rmvnorm_using_lower_chol(rng, this->mean + this->rho*(u - this->mean), this->lower_chol,
                                            std::sqrt(1.0 - this->rho*this->rho));
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

Parameters PCNProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters PCNProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const std::string &variable_in,
                                                                       const Particle &particle) const
{
  if (variable_in != this->variable)
    return Parameters();
  return this->simulate(rng, particle);
}

arma::mat PCNProposalKernel::specific_gradient_of_log(const std::string &variable_in,
                                                                            const Particle &proposed_particle,
                                                                            const Particle &old_particle)
{
  Rcpp::stop("ilike::pcn: the gradient of the proposal's log density is not available.");
}

arma::mat PCNProposalKernel::specific_subsample_gradient_of_log(const std::string &variable_in,
                                                                                      const Particle &proposed_particle,
                                                                                      const Particle &old_particle)
{
  Rcpp::stop("ilike::pcn: the gradient of the proposal's log density is not available.");
}

void PCNProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* PCNProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> PCNProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void PCNProposalKernel::set_index(Index* index_in)
{
}

void PCNProposalKernel::set_index_if_null(Index* index_in)
{
}

bool PCNProposalKernel::can_be_evaluated() const
{
  return true;
}

void PCNProposalKernel::set_data(Data* data_in)
{
}

}
