#include "custom_symmetric_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomSymmetricProposalKernel::CustomSymmetricProposalKernel()
  :SymmetricProposalKernel()
{
}

CustomSymmetricProposalKernel::~CustomSymmetricProposalKernel()
{
}

CustomSymmetricProposalKernel::CustomSymmetricProposalKernel(SimulateMCMCProposalPtr proposal_simulate_in)
  :SymmetricProposalKernel()
{
  //this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
}

CustomSymmetricProposalKernel::CustomSymmetricProposalKernel(const CustomSymmetricProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void CustomSymmetricProposalKernel::operator=(const CustomSymmetricProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomSymmetricProposalKernel::duplicate() const
{
  return( new CustomSymmetricProposalKernel(*this));
}

ProposalKernel* CustomSymmetricProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomSymmetricProposalKernel(*this));
}

SymmetricProposalKernel* CustomSymmetricProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new CustomSymmetricProposalKernel(*this));
}

void CustomSymmetricProposalKernel::make_copy(const CustomSymmetricProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
}

double CustomSymmetricProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle) const
{
  return this->proposal_evaluate(*proposed_particle.move_parameters,
                                 *old_particle.move_parameters,
                                 *this->proposal_parameters);
}

double CustomSymmetricProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters CustomSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 *particle.move_parameters,
                                 *this->proposal_parameters);
}


Parameters CustomSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}


Parameters CustomSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomSymmetricProposalKernel::subsample_simulate - not implemented.");
}

arma::mat CustomSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         Particle &proposed_particle,
                                                         Particle &old_particle)
{
  Rcpp::stop("CustomSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat CustomSymmetricProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                   Particle &proposed_particle,
                                                                   Particle &old_particle)
{
  Rcpp::stop("CustomSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomSymmetricProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}
