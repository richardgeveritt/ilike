#include "custom_distribution_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

CustomDistributionProposalKernel::CustomDistributionProposalKernel()
  :IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = NULL;
  //this->proposal_parameters = Parameters();
}

CustomDistributionProposalKernel::~CustomDistributionProposalKernel()
{
}

CustomDistributionProposalKernel::CustomDistributionProposalKernel(SimulateDistributionPtr proposal_simulate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = NULL;
  this->proposal_simulate = proposal_simulate_in;
}

CustomDistributionProposalKernel::CustomDistributionProposalKernel(SimulateDistributionPtr proposal_simulate_in,
                                                                   EvaluateLogDistributionPtr proposal_log_evaluate_in)
:IndependentProposalKernel()
{
  this->proposal_evaluate = proposal_log_evaluate_in;
  this->proposal_simulate = proposal_simulate_in;
}

CustomDistributionProposalKernel::CustomDistributionProposalKernel(const CustomDistributionProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void CustomDistributionProposalKernel::operator=(const CustomDistributionProposalKernel &another)
{
  if(this == &another)
    return;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CustomDistributionProposalKernel::duplicate() const
{
  return( new CustomDistributionProposalKernel(*this));
}

ProposalKernel* CustomDistributionProposalKernel::proposal_kernel_duplicate() const
{
  return( new CustomDistributionProposalKernel(*this));
}

IndependentProposalKernel* CustomDistributionProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new CustomDistributionProposalKernel(*this));
}

void CustomDistributionProposalKernel::make_copy(const CustomDistributionProposalKernel &another)
{
  this->proposal_evaluate = another.proposal_evaluate;
  this->proposal_simulate = another.proposal_simulate;
  //this->proposal_parameters = another.proposal_parameters;
}

double CustomDistributionProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  /*
  if (proposed_particle.is_empty())
    return this->proposal_evaluate(this->proposal_parameters);
  
  if (this->proposal_parameters.is_empty())
    return this->proposal_evaluate(proposed_particle);
  
  return this->proposal_evaluate(proposed_particle.merge(this->proposal_parameters));
  */
  return this->proposal_evaluate(proposed_particle);
}

/*
double CustomDistributionProposalKernel::evaluate_independent_kernel(Variables* proposed_particle,
                                                                     const Parameters &conditioned_on_parameters) const
{
  return this->proposal_evaluate(proposed_particle->move_parameters->merge(conditioned_on_parameters).merge(this->proposal_parameters));
}
*/

double CustomDistributionProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  /*
  // no difference since size of data set does not impact on proposal
  if (proposed_particle.is_empty())
  {
    this->proposal_evaluate(this->proposal_parameters);
  }
  
  if (this->proposal_parameters.is_empty())
  {
    this->proposal_evaluate(proposed_particle);
  }
  
  return this->proposal_evaluate(proposed_particle.merge(this->proposal_parameters));
  */
  
  return this->proposal_evaluate(proposed_particle);
}

/*
double CustomDistributionProposalKernel::subsample_evaluate_independent_kernel(Variables* proposed_particle,
                                                                               const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->proposal_evaluate(proposed_particle->move_parameters->merge(conditioned_on_parameters).merge(this->proposal_parameters));
}
*/

Parameters CustomDistributionProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

/*
Parameters CustomDistributionProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                 const Parameters &conditioned_on_parameters) const
{
  return this->proposal_simulate(rng,
                                 this->proposal_parameters);
}
*/

Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  // no difference since size of data set does not impact on proposal
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

/*
Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                            const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->proposal_simulate(rng,
                                 this->proposal_parameters);
}
*/

Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomDistributionProposalKernel::subsample_independent_simulate - not written yet.");
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

Parameters CustomDistributionProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                  const Parameters &conditioned_on_parameters) const
{
  /*
  if (this->proposal_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   conditioned_on_parameters);
  }
  
  if (conditioned_on_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters);
  }
  
  return this->proposal_simulate(rng,
                                 this->proposal_parameters.merge(conditioned_on_parameters));
  */
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                            const Parameters &conditioned_on_parameters) const
{
  /*
  // no difference since size of data set does not impact on proposal
  if (this->proposal_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   conditioned_on_parameters);
  }
  
  if (conditioned_on_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters);
  }
  
  return this->proposal_simulate(rng,
                                 this->proposal_parameters.merge(conditioned_on_parameters));
  */
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                            const std::string &variable,
                                                                            const Parameters &conditioned_on_parameters) const
{
  /*
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomDistributionProposalKernel::subsample_independent_simulate - not written yet.");
  if (this->proposal_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   conditioned_on_parameters);
  }
  
  if (conditioned_on_parameters.is_empty())
  {
    return this->proposal_simulate(rng,
                                   this->proposal_parameters);
  }
  
  return this->proposal_simulate(rng,
                                 this->proposal_parameters.merge(conditioned_on_parameters));
  */
  //return this->proposal_simulate(rng,
  //                               this->proposal_parameters);
  return this->proposal_simulate(rng);
}

/*
Parameters CustomDistributionProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                           const std::string &variable,
                                                                           const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  Rcpp::stop("CustomDistributionProposalKernel::subsample_independent_simulate - not written yet.");
  return this->proposal_simulate(rng,
                                 this->proposal_parameters);
}
*/

arma::mat CustomDistributionProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                        const Parameters &proposed_particle)
{
  Rcpp::stop("CustomDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}

/*
arma::mat CustomDistributionProposalKernel::independent_gradient_of_log(const std::string &variable,
                                      Variables* proposed_particle,
                                      const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}
*/

arma::mat CustomDistributionProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                  const Parameters &proposed_particle)
{
  Rcpp::stop("CustomDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}

void CustomDistributionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

/*
arma::mat CustomDistributionProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                Variables* proposed_particle,
                                                const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CustomDistributionProposalKernel::independent_gradient_of_log - not written yet.");
}
*/

GradientEstimatorOutput* CustomDistributionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CustomDistributionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void CustomDistributionProposalKernel::set_index(Index* index_in)
{
}
