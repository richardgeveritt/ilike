#include "truncated_normal_random_walk_proposal_kernel.h"
#include "truncated_normal_distribution_factor.h"

namespace ilike
{

TruncatedNormalRandomWalkProposalKernel::TruncatedNormalRandomWalkProposalKernel()
:ProposalKernel(), sd(1.0), lower(-arma::datum::inf), upper(arma::datum::inf)
{
}

TruncatedNormalRandomWalkProposalKernel::~TruncatedNormalRandomWalkProposalKernel()
{
}

TruncatedNormalRandomWalkProposalKernel::TruncatedNormalRandomWalkProposalKernel(const std::string &variable_in, double sd_in,
                                                                                 double lower_in, double upper_in)
:ProposalKernel(), variable(variable_in), sd(sd_in), lower(lower_in), upper(upper_in)
{
  if (!(sd_in > 0.0) || !(lower_in < upper_in))
    Rcpp::stop("ilike::tnorm_rw: the standard deviation must be positive and lower < upper.");
}

TruncatedNormalRandomWalkProposalKernel::TruncatedNormalRandomWalkProposalKernel(const TruncatedNormalRandomWalkProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void TruncatedNormalRandomWalkProposalKernel::operator=(const TruncatedNormalRandomWalkProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* TruncatedNormalRandomWalkProposalKernel::duplicate() const
{
  return new TruncatedNormalRandomWalkProposalKernel(*this);
}

ProposalKernel* TruncatedNormalRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return new TruncatedNormalRandomWalkProposalKernel(*this);
}

void TruncatedNormalRandomWalkProposalKernel::make_copy(const TruncatedNormalRandomWalkProposalKernel &another)
{
  this->variable = another.variable;
  this->sd = another.sd;
  this->lower = another.lower;
  this->upper = another.upper;
}

double TruncatedNormalRandomWalkProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                         const Particle &old_particle) const
{
  // log q(proposed | old): each element a truncated normal centred at its old value
  arma::colvec x = proposed_particle.get_transformed_parameters(this).get_colvec(this->variable);
  arma::colvec centre = old_particle.get_transformed_parameters(this).get_colvec(this->variable);
  if (x.n_elem != centre.n_elem)
    Rcpp::stop("ilike::tnorm_rw: the two points' blocks differ in size.");
  double output = 0.0;
  for (arma::uword k = 0; k < x.n_elem; ++k)
    output = output + log_truncated_normal(x[k], centre[k], this->sd, this->lower, this->upper);
  return output;
}

double TruncatedNormalRandomWalkProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                   const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters TruncatedNormalRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                             const Particle &particle) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x = current.get_colvec(this->variable);
  for (arma::uword k = 0; k < x.n_elem; ++k)
    x[k] = rtnorm(rng, x[k], this->sd, this->lower, this->upper);
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

// (C268) as simulate, written into the target (same draws, same values)
void TruncatedNormalRandomWalkProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                            const Particle &particle,
                                                            ProposalTarget &target) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x = current.get_colvec(this->variable);
  for (arma::uword k = 0; k < x.n_elem; ++k)
    x[k] = rtnorm(rng, x[k], this->sd, this->lower, this->upper);
  target.write(this->variable, x, current);
}

Parameters TruncatedNormalRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters TruncatedNormalRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                       const std::string &variable_in,
                                                                       const Particle &particle) const
{
  if (variable_in != this->variable)
    return Parameters();
  return this->simulate(rng, particle);
}

arma::mat TruncatedNormalRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable_in,
                                                                            const Particle &proposed_particle,
                                                                            const Particle &old_particle)
{
  Rcpp::stop("ilike::tnorm_rw: the gradient of the proposal's log density is not available.");
}

arma::mat TruncatedNormalRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable_in,
                                                                                      const Particle &proposed_particle,
                                                                                      const Particle &old_particle)
{
  Rcpp::stop("ilike::tnorm_rw: the gradient of the proposal's log density is not available.");
}

void TruncatedNormalRandomWalkProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* TruncatedNormalRandomWalkProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> TruncatedNormalRandomWalkProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void TruncatedNormalRandomWalkProposalKernel::set_index(Index* index_in)
{
}

void TruncatedNormalRandomWalkProposalKernel::set_index_if_null(Index* index_in)
{
}

bool TruncatedNormalRandomWalkProposalKernel::can_be_evaluated() const
{
  return true;
}

void TruncatedNormalRandomWalkProposalKernel::set_data(Data* data_in)
{
}

}
