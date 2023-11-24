#include <iterator>
#include "transformed_independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

TransformedIndependentProposalKernel::TransformedIndependentProposalKernel()
  :IndependentProposalKernel()
{
  this->proposal = NULL;
  this->distribution_on_transformed_space = true;
}

TransformedIndependentProposalKernel::~TransformedIndependentProposalKernel()
{
  if (this->proposal!=NULL)
    delete this->proposal;
}

// find mean and cov adaptively
TransformedIndependentProposalKernel::TransformedIndependentProposalKernel(IndependentProposalKernel* proposal_in,
                                                                           std::shared_ptr<Transform> transform_in,
                                                                           bool distribution_on_transformed_space_in)
  :IndependentProposalKernel(), transform(transform_in), distribution_on_transformed_space(distribution_on_transformed_space_in)
{
  this->proposal = proposal_in;
  
  if (this->transform==NULL)
  {
    Rcpp::stop("TransformedIndependentProposalKernel::TransformedIndependentProposalKernel - transform is NULL.");
  }
}

TransformedIndependentProposalKernel::TransformedIndependentProposalKernel(const TransformedIndependentProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void TransformedIndependentProposalKernel::operator=(const TransformedIndependentProposalKernel &another)
{
  if(this == &another)
    return;
  
  if (this->proposal!=NULL)
    delete this->proposal;

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* TransformedIndependentProposalKernel::duplicate() const
{
  return( new TransformedIndependentProposalKernel(*this));
}

ProposalKernel* TransformedIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new TransformedIndependentProposalKernel(*this));
}

IndependentProposalKernel* TransformedIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new TransformedIndependentProposalKernel(*this));
}

void TransformedIndependentProposalKernel::make_copy(const TransformedIndependentProposalKernel &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->independent_proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  this->transform = another.transform;
  this->distribution_on_transformed_space = another.distribution_on_transformed_space;
}

double TransformedIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  if (this->distribution_on_transformed_space==true)
  {
    return this->proposal->evaluate_independent_kernel(this->transform->inverse_transform(proposed_particle)) + this->transform->log_abs_inverse_jacobian_determinant(proposed_particle);
  }
  else
  {
    return this->proposal->evaluate_independent_kernel(proposed_particle);
  }
}

/*
double TransformedIndependentProposalKernel::evaluate_independent_kernel(Variables* proposed_particle,
                                                                      const Parameters &conditioned_on_parameters) const
{
  return this->evaluate_independent_kernel(proposed_particle);
}
*/

double TransformedIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  if (this->distribution_on_transformed_space==true)
  {
    return this->proposal->subsample_evaluate_independent_kernel(this->transform->inverse_transform(proposed_particle)) + this->transform->log_abs_inverse_jacobian_determinant(proposed_particle);
  }
  else
  {
    return this->proposal->subsample_evaluate_independent_kernel(proposed_particle);
  }
}

Parameters TransformedIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform->transform(this->proposal->independent_simulate(rng));
  else
  {
    Parameters output(this->proposal->independent_simulate(rng));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

Parameters TransformedIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                      const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
  {
    return this->transform->transform(this->proposal->independent_simulate(rng,
                                                                           conditioned_on_parameters));
  }
  else
  {
    Parameters output(this->proposal->independent_simulate(rng,
                                                           conditioned_on_parameters));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

Parameters TransformedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform->transform(this->proposal->subsample_independent_simulate(rng));
  else
  {
    Parameters output(this->proposal->independent_simulate(rng));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

Parameters TransformedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform->transform(this->proposal->subsample_independent_simulate(rng,
                                                                                     conditioned_on_parameters));
  else
  {
    Parameters output(this->proposal->independent_simulate(rng));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

Parameters TransformedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const std::string &variable) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform->transform(this->proposal->subsample_independent_simulate(rng,
                                                                                     variable));
  else
  {
    Parameters output(this->proposal->independent_simulate(rng));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

Parameters TransformedIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                const std::string &variable,
                                                                                const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform->transform(this->proposal->subsample_independent_simulate(rng,
                                                                                    variable,
                                                                                    conditioned_on_parameters));
  else
  {
    Parameters output(this->proposal->independent_simulate(rng));
    output.deep_overwrite_with_variables_in_argument(this->transform->transform(output));
    return output;
  }
}

arma::mat TransformedIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                            const Parameters &proposed_particle)
{
  Rcpp::stop("TransformedIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

arma::mat TransformedIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                      const Parameters &proposed_particle)
{
  Rcpp::stop("TransformedIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

void TransformedIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

GradientEstimatorOutput* TransformedIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> TransformedIndependentProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals = this->proposal->get_proposals();
  proposals.push_back(this);
  return proposals;
}

void TransformedIndependentProposalKernel::set_index(Index* index_in)
{
  this->proposal->set_index(index_in);
}

void TransformedIndependentProposalKernel::set_index_if_null(Index* index_in)
{
  this->proposal->set_index_if_null(index_in);
}
