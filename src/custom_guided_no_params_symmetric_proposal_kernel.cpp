#include "custom_guided_no_params_symmetric_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedNoParamsSymmetricProposalKernel::CustomGuidedNoParamsSymmetricProposalKernel()
  :SymmetricProposalKernel()
{
}

CustomGuidedNoParamsSymmetricProposalKernel::~CustomGuidedNoParamsSymmetricProposalKernel()
{
}

CustomGuidedNoParamsSymmetricProposalKernel::CustomGuidedNoParamsSymmetricProposalKernel(SimulateGuidedNoParamsMCMCProposalPtr proposal_simulate_in,
                                                                       const Data* data_in)
  :SymmetricProposalKernel()
{
  //this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedNoParamsSymmetricProposalKernel::CustomGuidedNoParamsSymmetricProposalKernel(const CustomGuidedNoParamsSymmetricProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedNoParamsSymmetricProposalKernel::operator=(const CustomGuidedNoParamsSymmetricProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedNoParamsSymmetricProposalKernel::duplicate() const
{
  return( new CustomGuidedNoParamsSymmetricProposalKernel(*this));
}

ProposalKernel* CustomGuidedNoParamsSymmetricProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedNoParamsSymmetricProposalKernel(*this));
}

SymmetricProposalKernel* CustomGuidedNoParamsSymmetricProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new CustomGuidedNoParamsSymmetricProposalKernel(*this));
}

void CustomGuidedNoParamsSymmetricProposalKernel::make_copy(const CustomGuidedNoParamsSymmetricProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->data = another.data;
}

double CustomGuidedNoParamsSymmetricProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                             const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this),
                                 *this->data);
}

/*
double CustomGuidedNoParamsSymmetricProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.move_parameters->merge(conditioned_on_parameters),
                                 old_particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomGuidedNoParamsSymmetricProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                                       const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle.get_transformed_parameters(this),
                                        old_particle.get_transformed_parameters(this));
}

/*
double CustomGuidedNoParamsSymmetricProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle, conditioned_on_parameters);
}
*/

Parameters CustomGuidedNoParamsSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                                                 const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this),
                                 *this->data);
}

/*
Parameters CustomGuidedNoParamsSymmetricProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

Parameters CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                           const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
Parameters CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle, conditioned_on_parameters);
}
*/

Parameters CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable,
                                                                           const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate - not implemented.");
}

/*
Parameters CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedNoParamsSymmetricProposalKernel::subsample_simulate - not implemented.");
}
*/

arma::mat CustomGuidedNoParamsSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                                const Particle &proposed_particle,
                                                                                const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat CustomGuidedNoParamsSymmetricProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat CustomGuidedNoParamsSymmetricProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                          const Particle &proposed_particle,
                                                                                          const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedNoParamsSymmetricProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomGuidedNoParamsSymmetricProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{

}

GradientEstimatorOutput* CustomGuidedNoParamsSymmetricProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomGuidedNoParamsSymmetricProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomGuidedNoParamsSymmetricProposalKernel::set_index(Index* index_in)
{
}
