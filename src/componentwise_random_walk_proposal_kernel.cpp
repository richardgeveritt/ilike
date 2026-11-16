#include <cmath>
#include "componentwise_random_walk_proposal_kernel.h"
#include "distributions.h"

namespace ilike
{

ComponentwiseRandomWalkProposalKernel::ComponentwiseRandomWalkProposalKernel()
:GaussianRandomWalkProposalKernel()
{
}

ComponentwiseRandomWalkProposalKernel::~ComponentwiseRandomWalkProposalKernel()
{
}

ComponentwiseRandomWalkProposalKernel::ComponentwiseRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                             const arma::mat &covariance_in,
                                                                             double scale_in)
:GaussianRandomWalkProposalKernel(variable_name_in, covariance_in, scale_in)
{
  this->variable = variable_name_in;
  this->component_scales = arma::colvec(covariance_in.n_rows, arma::fill::value(scale_in));
}

ComponentwiseRandomWalkProposalKernel::ComponentwiseRandomWalkProposalKernel(const ComponentwiseRandomWalkProposalKernel &another)
:GaussianRandomWalkProposalKernel(another)
{
  this->make_copy(another);
}

void ComponentwiseRandomWalkProposalKernel::operator=(const ComponentwiseRandomWalkProposalKernel &another)
{
  if(this == &another)
    return;
  GaussianRandomWalkProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* ComponentwiseRandomWalkProposalKernel::duplicate() const
{
  return( new ComponentwiseRandomWalkProposalKernel(*this));
}

ProposalKernel* ComponentwiseRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return( new ComponentwiseRandomWalkProposalKernel(*this));
}

SymmetricProposalKernel* ComponentwiseRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new ComponentwiseRandomWalkProposalKernel(*this));
}

void ComponentwiseRandomWalkProposalKernel::make_copy(const ComponentwiseRandomWalkProposalKernel &another)
{
  this->variable = another.variable;
  this->component_scales = another.component_scales;
}

arma::colvec* ComponentwiseRandomWalkProposalKernel::get_component_scales(const std::string &variable_in)
{
  if (variable_in!=this->variable)
    return NULL;
  return &this->component_scales;
}

arma::colvec ComponentwiseRandomWalkProposalKernel::move_one_component(RandomNumberGenerator &rng,
                                                                      const Particle &particle) const
{
  arma::colvec point = particle.get_transformed_parameters(this).get_colvec(this->variable);
  size_t d = point.n_elem;
  if (d!=this->component_scales.n_elem)
    Rcpp::stop("ilike::componentwise_rw: the covariance does not match the dimension of " + this->variable + ".");
  // k uniform on 0, ..., d-1, then one normal increment for component k
  size_t k = size_t(std::floor(runif(rng)*double(d)));
  if (k>=d)
    k = d-1;
  const GaussianProposalInfo &info = this->proposal_info.find(this->variable)->second;
  point[k] = point[k] + std::sqrt(this->component_scales[k]*info.get_covariance()(k,k))*rnorm(rng);
  return point;
}

Parameters ComponentwiseRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                           const Particle &particle) const
{
  Parameters output;
  output[this->variable] = this->move_one_component(rng, particle);
  return output;
}

void ComponentwiseRandomWalkProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                          const Particle &particle,
                                                          ProposalTarget &target) const
{
  target.matrix(this->variable) = this->move_one_component(rng, particle);
}

Parameters ComponentwiseRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                     const Particle &particle) const
{
  return this->simulate(rng, particle);
}

double ComponentwiseRandomWalkProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                       const Particle &old_particle) const
{
  Rcpp::stop("ilike::componentwise_rw: its density is singular (one component moves); use it in an m_proposal, where it cancels.");
}

double ComponentwiseRandomWalkProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                                       const Parameters &old_parameters) const
{
  Rcpp::stop("ilike::componentwise_rw: its density is singular (one component moves); use it in an m_proposal, where it cancels.");
}

double ComponentwiseRandomWalkProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                 const Particle &old_particle) const
{
  Rcpp::stop("ilike::componentwise_rw: its density is singular (one component moves); use it in an m_proposal, where it cancels.");
}

}
