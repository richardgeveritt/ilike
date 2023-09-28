#include "custom_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomProposalKernel::CustomProposalKernel()
  :ProposalKernel()
{
}

CustomProposalKernel::~CustomProposalKernel()
{
}

CustomProposalKernel::CustomProposalKernel(SimulateMCMCProposalPtr proposal_simulate_in,
                                           EvaluateLogMCMCProposalPtr proposal_evaluate_in)
  :ProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
}

CustomProposalKernel::CustomProposalKernel(const CustomProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void CustomProposalKernel::operator=(const CustomProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomProposalKernel::duplicate() const
{
  return( new CustomProposalKernel(*this));
}

ProposalKernel* CustomProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomProposalKernel(*this));
}

void CustomProposalKernel::make_copy(const CustomProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
}

double CustomProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                      const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this),
                                 *this->proposal_parameters);
}

double CustomProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters CustomProposalKernel::simulate(RandomNumberGenerator &rng,
                                          const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this),
                                 *this->proposal_parameters);
}


Parameters CustomProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}


Parameters CustomProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomProposalKernel::subsample_simulate - not implemented.");
}

arma::mat CustomProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         const Particle &proposed_particle,
                                                         const Particle &old_particle)
{
  Rcpp::stop("CustomProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat CustomProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                   const Particle &proposed_particle,
                                                                   const Particle &old_particle)
{
  Rcpp::stop("CustomProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}

GradientEstimatorOutput* CustomProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomProposalKernel::set_index(Index* index_in)
{
}

void CustomProposalKernel::set_index_if_null(Index* index_in)
{
}
