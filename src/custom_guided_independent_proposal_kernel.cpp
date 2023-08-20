#include "custom_guided_independent_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedIndependentProposalKernel::CustomGuidedIndependentProposalKernel()
  :IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = NULL;
  this->proposal_parameters = NULL;
  this->data = NULL;
}

CustomGuidedIndependentProposalKernel::~CustomGuidedIndependentProposalKernel()
{
}

CustomGuidedIndependentProposalKernel::CustomGuidedIndependentProposalKernel(SimulateGuidedIndependentProposalPtr proposal_simulate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
  this->proposal_parameters = NULL;
  this->data = NULL;
}

/*
CustomGuidedIndependentProposalKernel::CustomGuidedIndependentProposalKernel(SimulateGuidedIndependentProposalPtr proposal_simulate_in,
                                                                 EvaluateLogGuidedIndependentProposalPtr proposal_evaluate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->proposal_parameters = NULL;
  this->data = NULL;
}
*/

CustomGuidedIndependentProposalKernel::CustomGuidedIndependentProposalKernel(SimulateGuidedIndependentProposalPtr proposal_simulate_in,
                                                                 EvaluateLogGuidedIndependentProposalPtr proposal_evaluate_in,
                                                                             const Data* data_in)
  :IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->proposal_parameters = NULL;
  this->data = data_in;
}

CustomGuidedIndependentProposalKernel::CustomGuidedIndependentProposalKernel(const CustomGuidedIndependentProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedIndependentProposalKernel::operator=(const CustomGuidedIndependentProposalKernel &another)
{
  if(this == &another)
    return;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedIndependentProposalKernel::duplicate() const
{
  return( new CustomGuidedIndependentProposalKernel(*this));
}

ProposalKernel* CustomGuidedIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedIndependentProposalKernel(*this));
}

IndependentProposalKernel* CustomGuidedIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new CustomGuidedIndependentProposalKernel(*this));
}

void CustomGuidedIndependentProposalKernel::make_copy(const CustomGuidedIndependentProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
  this->data = another.data;
}

void CustomGuidedIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}

double CustomGuidedIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_evaluate(proposed_particle,
                                   *this->proposal_parameters,
                                   *this->data);
  }
  else
  {
    stop("CustomGuidedIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

/*
double CustomGuidedIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle,
                                                                    const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomGuidedIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  stop("CustomGuidedIndependentProposalKernel::subsample_evaluate_independent_kernel - not yet written.");
}

Parameters CustomGuidedIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   *this->proposal_parameters,
                                   *this->data);
  }
  else
  {
    stop("CustomGuidedIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomGuidedIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                 const Parameters &conditioned_on_parameters) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters->merge(conditioned_on_parameters),
                                   *this->data);
  }
  else
  {
    stop("CustomGuidedIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomGuidedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  stop("CustomGuidedIndependentProposalKernel::subsample_independent_simulate - not yet written.");
  // no difference since size of data set does not impact on proposal
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
}

Parameters CustomGuidedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const Parameters &conditioned_on_parameters) const
{
  stop("CustomGuidedIndependentProposalKernel::subsample_independent_simulate - not yet written.");
  // no difference since size of data set does not impact on proposal
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters.merge(conditioned_on_parameters));
}

Parameters CustomGuidedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedIndependentProposalKernel::subsample_independent_simulate - not written yet.");
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
}

Parameters CustomGuidedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable,
                                                                           const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedIndependentProposalKernel::subsample_independent_simulate - not written yet.");
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters.merge(conditioned_on_parameters));
}

arma::mat CustomGuidedIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                       const Parameters &proposed_particle)
{
  Rcpp::stop("CustomGuidedIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

/*
arma::mat CustomGuidedIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                       Variables* proposed_particle,
                                                                       const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}
*/

arma::mat CustomGuidedIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                 const Parameters &proposed_particle)
{
  Rcpp::stop("CustomGuidedIndependentProposalKernel::subsample_independent_gradient_of_log - not written yet.");
}

/*
arma::mat CustomGuidedIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                             Variables* proposed_particle,
                                             const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedIndependentProposalKernel::subsample_independent_gradient_of_log - not written yet.");
}
*/

GradientEstimatorOutput* CustomGuidedIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<ProposalKernel*> CustomGuidedIndependentProposalKernel::get_proposals()
{
  std::vector<ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomGuidedIndependentProposalKernel::set_index(Index* index_in)
{
}
