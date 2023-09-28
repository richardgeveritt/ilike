#include "custom_no_params_symmetric_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomNoParamsSymmetricProposalKernel::CustomNoParamsSymmetricProposalKernel()
  :SymmetricProposalKernel()
{
}

CustomNoParamsSymmetricProposalKernel::~CustomNoParamsSymmetricProposalKernel()
{
}

CustomNoParamsSymmetricProposalKernel::CustomNoParamsSymmetricProposalKernel(SimulateNoParamsMCMCProposalPtr proposal_simulate_in)
  :SymmetricProposalKernel()
{
  //this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
}

CustomNoParamsSymmetricProposalKernel::CustomNoParamsSymmetricProposalKernel(const CustomNoParamsSymmetricProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void CustomNoParamsSymmetricProposalKernel::operator=(const CustomNoParamsSymmetricProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomNoParamsSymmetricProposalKernel::duplicate() const
{
  return( new CustomNoParamsSymmetricProposalKernel(*this));
}

ProposalKernel* CustomNoParamsSymmetricProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomNoParamsSymmetricProposalKernel(*this));
}

SymmetricProposalKernel* CustomNoParamsSymmetricProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new CustomNoParamsSymmetricProposalKernel(*this));
}

void CustomNoParamsSymmetricProposalKernel::make_copy(const CustomNoParamsSymmetricProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
}

double CustomNoParamsSymmetricProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                       const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this));
}

/*
double CustomNoParamsSymmetricProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.move_parameters->merge(conditioned_on_parameters),
                                 old_particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomNoParamsSymmetricProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                 const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
double CustomNoParamsSymmetricProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle, conditioned_on_parameters);
}
*/

Parameters CustomNoParamsSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                                           const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this));
}

/*
Parameters CustomNoParamsSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

Parameters CustomNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                     const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
Parameters CustomNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle, conditioned_on_parameters);
}
*/

Parameters CustomNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                     const std::string &variable,
                                                                     const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomNoParamsSymmetricProposalKernel::subsample_simulate - not implemented.");
}

/*
Parameters CustomNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomNoParamsSymmetricProposalKernel::subsample_simulate - not implemented.");
}
*/

arma::mat CustomNoParamsSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                          const Particle &proposed_particle,
                                                                          const Particle &old_particle)
{
  Rcpp::stop("CustomNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat CustomNoParamsSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat CustomNoParamsSymmetricProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                    const Particle &proposed_particle,
                                                                                    const Particle &old_particle)
{
  Rcpp::stop("CustomNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomNoParamsSymmetricProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{

}

GradientEstimatorOutput* CustomNoParamsSymmetricProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomNoParamsSymmetricProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomNoParamsSymmetricProposalKernel::set_index(Index* index_in)
{
}

void CustomNoParamsSymmetricProposalKernel::set_index_if_null(Index* index_in)
{
}
