#include "custom_no_params_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomNoParamsProposalKernel::CustomNoParamsProposalKernel()
  :ProposalKernel()
{
}

CustomNoParamsProposalKernel::~CustomNoParamsProposalKernel()
{
}

CustomNoParamsProposalKernel::CustomNoParamsProposalKernel(SimulateNoParamsMCMCProposalPtr proposal_simulate_in,
                                                           EvaluateLogNoParamsMCMCProposalPtr proposal_evaluate_in)
  :ProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
}

CustomNoParamsProposalKernel::CustomNoParamsProposalKernel(const CustomNoParamsProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void CustomNoParamsProposalKernel::operator=(const CustomNoParamsProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomNoParamsProposalKernel::duplicate() const
{
  return( new CustomNoParamsProposalKernel(*this));
}

ProposalKernel* CustomNoParamsProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomNoParamsProposalKernel(*this));
}

void CustomNoParamsProposalKernel::make_copy(const CustomNoParamsProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
}

double CustomNoParamsProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                              const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this));
}

/*
double CustomNoParamsProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.move_parameters->merge(conditioned_on_parameters),
                                 old_particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomNoParamsProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                        const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
double CustomNoParamsProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle, conditioned_on_parameters);
}
*/

Parameters CustomNoParamsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                  const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this));
}

/*
Parameters CustomNoParamsProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

Parameters CustomNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
Parameters CustomNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle, conditioned_on_parameters);
}
*/

Parameters CustomNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            const std::string &variable,
                                                            const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomNoParamsProposalKernel::subsample_simulate - not implemented.");
}

/*
Parameters CustomNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomNoParamsProposalKernel::subsample_simulate - not implemented.");
}
*/

arma::mat CustomNoParamsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                 const Particle &proposed_particle,
                                                                 const Particle &old_particle)
{
  Rcpp::stop("CustomNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat CustomNoParamsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat CustomNoParamsProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                           const Particle &proposed_particle,
                                                                           const Particle &old_particle)
{
  Rcpp::stop("CustomNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomNoParamsProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{

}

GradientEstimatorOutput* CustomNoParamsProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomNoParamsProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomNoParamsProposalKernel::set_index(Index* index_in)
{
}

void CustomNoParamsProposalKernel::set_index_if_null(Index* index_in)
{
}
