#include "custom_guided_symmetric_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedSymmetricProposalKernel::CustomGuidedSymmetricProposalKernel()
  :SymmetricProposalKernel()
{
}

CustomGuidedSymmetricProposalKernel::~CustomGuidedSymmetricProposalKernel()
{
}

CustomGuidedSymmetricProposalKernel::CustomGuidedSymmetricProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                                                       const Data* data_in)
  :SymmetricProposalKernel()
{
  //this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedSymmetricProposalKernel::CustomGuidedSymmetricProposalKernel(const CustomGuidedSymmetricProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedSymmetricProposalKernel::operator=(const CustomGuidedSymmetricProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedSymmetricProposalKernel::duplicate() const
{
  return( new CustomGuidedSymmetricProposalKernel(*this));
}

ProposalKernel* CustomGuidedSymmetricProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedSymmetricProposalKernel(*this));
}

SymmetricProposalKernel* CustomGuidedSymmetricProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new CustomGuidedSymmetricProposalKernel(*this));
}

void CustomGuidedSymmetricProposalKernel::make_copy(const CustomGuidedSymmetricProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
  this->data = another.data;
}

double CustomGuidedSymmetricProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                            Particle &old_particle) const
{
  return this->proposal_evaluate(*proposed_particle.move_parameters,
                                 *old_particle.move_parameters,
                                 *this->proposal_parameters,
                                 *this->data);
}

double CustomGuidedSymmetricProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                      Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters CustomGuidedSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                                Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 *particle.move_parameters,
                                 *this->proposal_parameters,
                                 *this->data);
}


Parameters CustomGuidedSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                          Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}


Parameters CustomGuidedSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                          const std::string &variable,
                                                          Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedSymmetricProposalKernel::subsample_simulate - not implemented.");
}

arma::mat CustomGuidedSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                               Particle &proposed_particle,
                                                               Particle &old_particle)
{
  Rcpp::stop("CustomGuidedSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat CustomGuidedSymmetricProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                         Particle &proposed_particle,
                                                                         Particle &old_particle)
{
  Rcpp::stop("CustomGuidedSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomGuidedSymmetricProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}
