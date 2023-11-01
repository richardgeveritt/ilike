#include "custom_guided_no_params_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedNoParamsProposalKernel::CustomGuidedNoParamsProposalKernel()
  :ProposalKernel()
{
}

CustomGuidedNoParamsProposalKernel::~CustomGuidedNoParamsProposalKernel()
{
}

CustomGuidedNoParamsProposalKernel::CustomGuidedNoParamsProposalKernel(SimulateGuidedNoParamsMCMCProposalPtr proposal_simulate_in,
                                                                       const Data* data_in)
:ProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedNoParamsProposalKernel::CustomGuidedNoParamsProposalKernel(SimulateGuidedNoParamsMCMCProposalPtr proposal_simulate_in,
                                                                       EvaluateLogGuidedNoParamsMCMCProposalPtr proposal_evaluate_in,
                                                                       const Data* data_in)
  :ProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedNoParamsProposalKernel::CustomGuidedNoParamsProposalKernel(const CustomGuidedNoParamsProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedNoParamsProposalKernel::operator=(const CustomGuidedNoParamsProposalKernel &another)
{
  if(this == &another)
    return;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedNoParamsProposalKernel::duplicate() const
{
  return( new CustomGuidedNoParamsProposalKernel(*this));
}

ProposalKernel* CustomGuidedNoParamsProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedNoParamsProposalKernel(*this));
}

void CustomGuidedNoParamsProposalKernel::make_copy(const CustomGuidedNoParamsProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->data = another.data;
}

double CustomGuidedNoParamsProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                    const Particle &old_particle) const
{
  return this->proposal_evaluate(proposed_particle.get_transformed_parameters(this),
                                 old_particle.get_transformed_parameters(this),
                                 *this->data);
}

/*
double CustomGuidedNoParamsProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                      Particle &old_particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.move_parameters->merge(conditioned_on_parameters),
                                 old_particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomGuidedNoParamsProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                              const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
double CustomGuidedNoParamsProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                Particle &old_particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle, conditioned_on_parameters);
}
*/

Parameters CustomGuidedNoParamsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                        const Particle &particle) const
{
  return this->proposal_simulate(rng,
                                 particle.get_transformed_parameters(this),
                                 *this->data);
}

/*
Parameters CustomGuidedNoParamsProposalKernel::simulate(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 particle.move_parameters->merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

Parameters CustomGuidedNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                  const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
Parameters CustomGuidedNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle, conditioned_on_parameters);
}
*/

Parameters CustomGuidedNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                  const std::string &variable,
                                                                  const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedNoParamsProposalKernel::subsample_simulate - not implemented.");
}

/*
Parameters CustomGuidedNoParamsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                    const std::string &variable,
                                                    Particle &particle,
                                                    const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedNoParamsProposalKernel::subsample_simulate - not implemented.");
}
*/

arma::mat CustomGuidedNoParamsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                       const Particle &proposed_particle,
                                                                       const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat CustomGuidedNoParamsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat CustomGuidedNoParamsProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                                 const Particle &proposed_particle,
                                                                                 const Particle &old_particle)
{
  Rcpp::stop("CustomGuidedNoParamsProposalKernel::specific_gradient_of_log - not written yet.");
}

void CustomGuidedNoParamsProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{

}

GradientEstimatorOutput* CustomGuidedNoParamsProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomGuidedNoParamsProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomGuidedNoParamsProposalKernel::set_index(Index* index_in)
{
}

void CustomGuidedNoParamsProposalKernel::set_index_if_null(Index* index_in)
{
}
