#include "custom_guided_distribution_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomGuidedDistributionProposalKernel::CustomGuidedDistributionProposalKernel()
  :IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = NULL;
  this->data = NULL;
}

CustomGuidedDistributionProposalKernel::~CustomGuidedDistributionProposalKernel()
{
}

CustomGuidedDistributionProposalKernel::CustomGuidedDistributionProposalKernel(SimulateGuidedDistributionPtr proposal_simulate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
  this->data = NULL;
}

/*
CustomGuidedDistributionProposalKernel::CustomGuidedDistributionProposalKernel(SimulateGuidedDistributionPtr proposal_simulate_in,
                                                                 EvaluateLogGuidedDistributionPtr proposal_evaluate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = NULL;
}
*/

CustomGuidedDistributionProposalKernel::CustomGuidedDistributionProposalKernel(SimulateGuidedDistributionPtr proposal_simulate_in,
                                                                 EvaluateLogGuidedDistributionPtr proposal_evaluate_in,
                                                                             const Data* data_in)
  :IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->data = data_in;
}

CustomGuidedDistributionProposalKernel::CustomGuidedDistributionProposalKernel(const CustomGuidedDistributionProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void CustomGuidedDistributionProposalKernel::operator=(const CustomGuidedDistributionProposalKernel &another)
{
  if(this == &another)
    return;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomGuidedDistributionProposalKernel::duplicate() const
{
  return( new CustomGuidedDistributionProposalKernel(*this));
}

ProposalKernel* CustomGuidedDistributionProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomGuidedDistributionProposalKernel(*this));
}

IndependentProposalKernel* CustomGuidedDistributionProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new CustomGuidedDistributionProposalKernel(*this));
}

void CustomGuidedDistributionProposalKernel::make_copy(const CustomGuidedDistributionProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->data = another.data;
}

double CustomGuidedDistributionProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  return this->proposal_evaluate(proposed_particle,
                                 *this->data);
}

/*
double CustomGuidedDistributionProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle,
                                                                    const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle.merge(conditioned_on_parameters),
                                 this->proposal_parameters);
}
*/

double CustomGuidedDistributionProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  stop("CustomGuidedDistributionProposalKernel::subsample_evaluate_independent_kernel - not yet written.");
}

Parameters CustomGuidedDistributionProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  return this->proposal_simulate(rng,
                                 *this->data);
}

Parameters CustomGuidedDistributionProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                        const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 *this->data);
}

Parameters CustomGuidedDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  stop("CustomGuidedDistributionProposalKernel::subsample_independent_simulate - not yet written.");
  // no difference since size of data set does not impact on proposal
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
}

Parameters CustomGuidedDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                  const Parameters &conditioned_on_parameters) const
{
  stop("CustomGuidedDistributionProposalKernel::subsample_independent_simulate - not yet written.");
  // no difference since size of data set does not impact on proposal
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters.merge(conditioned_on_parameters));
}

Parameters CustomGuidedDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                  const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedDistributionProposalKernel::subsample_independent_simulate - not written yet.");
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
}

Parameters CustomGuidedDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable,
                                                                           const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomGuidedDistributionProposalKernel::subsample_independent_simulate - not written yet.");
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters.merge(conditioned_on_parameters));
}

arma::mat CustomGuidedDistributionProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                       const Parameters &proposed_particle)
{
  Rcpp::stop("CustomGuidedDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}

/*
arma::mat CustomGuidedDistributionProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                       Variables* proposed_particle,
                                                                       const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}
*/

arma::mat CustomGuidedDistributionProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                 const Parameters &proposed_particle)
{
  Rcpp::stop("CustomGuidedDistributionProposalKernel::subsample_independent_gradient_of_log - not written yet.");
}

void CustomGuidedDistributionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

/*
arma::mat CustomGuidedDistributionProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                             Variables* proposed_particle,
                                             const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomGuidedDistributionProposalKernel::subsample_independent_gradient_of_log - not written yet.");
}
*/

GradientEstimatorOutput* CustomGuidedDistributionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomGuidedDistributionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomGuidedDistributionProposalKernel::set_index(Index* index_in)
{
}
