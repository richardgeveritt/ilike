#include "custom_guided_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedProposalKernel::CustomGuidedProposalKernel()
  :ProposalKernel()
{
}

CustomGuidedProposalKernel::~CustomGuidedProposalKernel()
{
}

CustomGuidedProposalKernel::CustomGuidedProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                                                       const Data* data_in)
:ProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedProposalKernel::CustomGuidedProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                                                       EvaluateLogGuidedMCMCProposalPtr proposal_evaluate_in,
                                                       const Data* data_in)
  :ProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedProposalKernel::CustomGuidedProposalKernel(const CustomGuidedProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedProposalKernel::operator=(const CustomGuidedProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedProposalKernel::duplicate() const
{
  return( new CustomGuidedProposalKernel(*this));
}

ProposalKernel* CustomGuidedProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedProposalKernel(*this));
}

void CustomGuidedProposalKernel::make_copy(const CustomGuidedProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
  this->data = another.data;
}

double CustomGuidedProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                            const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this),
                                 *this->proposal_parameters,
                                 *this->data);
}

double CustomGuidedProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                      const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters CustomGuidedProposalKernel::simulate(RandomNumberGenerator &rng,
                                                const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this),
                                 *this->proposal_parameters,
                                 *this->data);
}


Parameters CustomGuidedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                          const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}


Parameters CustomGuidedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                          const std::string &variable,
                                                          const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedProposalKernel::subsample_simulate - not implemented.");
}

arma::mat CustomGuidedProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                               const Particle &proposed_particle,
                                                               const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat CustomGuidedProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                         const Particle &proposed_particle,
                                                                         const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomGuidedProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}

GradientEstimatorOutput* CustomGuidedProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomGuidedProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomGuidedProposalKernel::set_index(Index* index_in)
{
}

void CustomGuidedProposalKernel::set_index_if_null(Index* index_in)
{
}
