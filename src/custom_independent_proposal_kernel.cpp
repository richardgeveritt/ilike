#include "custom_independent_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

namespace ilike
{
CustomIndependentProposalKernel::CustomIndependentProposalKernel()
:IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = NULL;
  this->proposal_parameters = NULL;
}

CustomIndependentProposalKernel::~CustomIndependentProposalKernel()
{
}

CustomIndependentProposalKernel::CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
  this->proposal_parameters = NULL;
}

CustomIndependentProposalKernel::CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
                                                                 EvaluateLogIndependentProposalPtr proposal_evaluate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
  this->proposal_parameters = NULL;
}

/*
 CustomIndependentProposalKernel::CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
 EvaluateLogIndependentProposalPtr proposal_evaluate_in,
 Parameters* proposal_parameters_in)
 :IndependentProposalKernel()
 {
 this->proposal_evaluate = proposal_evaluate_in;
 this->proposal_simulate = proposal_simulate_in;
 this->proposal_parameters = proposal_parameters_in;
 }
 */

CustomIndependentProposalKernel::CustomIndependentProposalKernel(const CustomIndependentProposalKernel &another)
:IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void CustomIndependentProposalKernel::operator=(const CustomIndependentProposalKernel &another)
{
  if(this == &another)
    return;
  
  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomIndependentProposalKernel::duplicate() const
{
  return( new CustomIndependentProposalKernel(*this));
}

ProposalKernel* CustomIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomIndependentProposalKernel(*this));
}

IndependentProposalKernel* CustomIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new CustomIndependentProposalKernel(*this));
}

void CustomIndependentProposalKernel::make_copy(const CustomIndependentProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  this->proposal_parameters = another.proposal_parameters;
}

void CustomIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal_parameters = proposal_parameters_in;
}

double CustomIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_evaluate(proposed_particle,
                                   *this->proposal_parameters);
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

/*
 double CustomIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle,
 const Parameters &conditioned_on_parameters) const
 {
 return this->proposal_evaluate(proposed_particle.merge(conditioned_on_parameters),
 this->proposal_parameters);
 }
 */

double CustomIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_evaluate(proposed_particle,
                                   *this->proposal_parameters);
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   *this->proposal_parameters);
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                 const Parameters &conditioned_on_parameters) const
{
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters->merge(conditioned_on_parameters));
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  // no difference since size of data set does not impact on proposal
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   *this->proposal_parameters);
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  if (this->proposal_parameters!=NULL)
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters->merge(conditioned_on_parameters));
  }
  else
  {
    stop("CustomIndependentProposalKernel::evaluate_independent_kernel - proposal parameters not set.");
  }
}

Parameters CustomIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomIndependentProposalKernel::subsample_independent_simulate - not written yet.");
  return this->proposal_simulate(rng,
                                 *this->proposal_parameters);
}

Parameters CustomIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable,
                                                                           const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomIndependentProposalKernel::subsample_independent_simulate - not written yet.");
  return this->proposal_simulate(rng,
                                 this->proposal_parameters->merge(conditioned_on_parameters));
}

arma::mat CustomIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                       const Parameters &proposed_particle)
{
  Rcpp::stop("CustomIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

/*
 arma::mat CustomIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
 Variables* proposed_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("CustomIndependentProposalKernel::independent_gradient_of_log - not written yet.");
 }
 */

arma::mat CustomIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                 const Parameters &proposed_particle)
{
  Rcpp::stop("CustomIndependentProposalKernel::subsample_independent_gradient_of_log - not written yet.");
}

/*
 arma::mat CustomIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
 Variables* proposed_particle,
 const Parameters &conditioned_on_parameters)
 {
 Rcpp::stop("CustomIndependentProposalKernel::subsample_independent_gradient_of_log - not written yet.");
 }
 */

GradientEstimatorOutput* CustomIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomIndependentProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomIndependentProposalKernel::set_index(Index* index_in)
{
}

void CustomIndependentProposalKernel::set_index_if_null(Index* index_in)
{
}

bool CustomIndependentProposalKernel::can_be_evaluated() const
{
  if (this->proposal_evaluate==NULL)
    return false;
  else
    return true;
}

void CustomIndependentProposalKernel::set_data(Data* data_in)
{
  
}
}
