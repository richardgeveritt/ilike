#include "transformed_proposal_kernel.h"
#include "mcmc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "smc_adaptor.h"

TransformedProposalKernel::TransformedProposalKernel()
  :ProposalKernel()
{
  this->proposal = NULL;
  this->distribution_on_transformed_space = true;
}

TransformedProposalKernel::~TransformedProposalKernel()
{
  if (this->proposal!=NULL)
    delete this->proposal;
}

TransformedProposalKernel::TransformedProposalKernel(ProposalKernel* proposal_in,
                                                     const Transform &transform_in,
                                                     bool distribution_on_transformed_space_in)
:ProposalKernel(), transform(transform_in), distribution_on_transformed_space(distribution_on_transformed_space_in)
{
  this->proposal = proposal_in;
}

TransformedProposalKernel::TransformedProposalKernel(const TransformedProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void TransformedProposalKernel::operator=(const TransformedProposalKernel &another)
{
  if(this == &another)
    return;
  
  if (this->proposal!=NULL)
    delete this->proposal;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* TransformedProposalKernel::duplicate() const
{
  return( new TransformedProposalKernel(*this));
}

ProposalKernel* TransformedProposalKernel::proposal_kernel_duplicate() const
{
  return( new TransformedProposalKernel(*this));
}

void TransformedProposalKernel::make_copy(const TransformedProposalKernel &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  this->transform = another.transform;
  this->distribution_on_transformed_space = another.distribution_on_transformed_space;
}

double TransformedProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                           const Particle &old_particle) const
{
  Rcpp::stop("TransformedProposalKernel::specific_evaluate_kernel - option for distribution_on_transformed_space not written yet.");
  
  if (this->distribution_on_transformed_space==true)
  {
    //proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform.inverse_transform(proposed_particle.parameters));
    //return this->proposal->evaluate_kernel(proposed_particle,old_particle) + this->transform.log_abs_inverse_jacobian_determinant(proposed_particle.parameters);
  }
  else
  {
    return this->proposal->evaluate_kernel(proposed_particle,old_particle);
  }
}

/*
double TransformedProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
  {
    proposed_particle.parameters.add_parameters(this->transform.inverse_transform(proposed_particle.parameters));
    return this->proposal->evaluate_kernel(proposed_particle,
                                           old_particle,
                                           conditioned_on_parameters) + this->transform.log_abs_inverse_jacobian_determinant(proposed_particle.parameters);
  }
  else
  {
    return this->proposal->evaluate_kernel(proposed_particle,old_particle,conditioned_on_parameters);
  }
}
*/

double TransformedProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                     const Particle &old_particle) const
{
  Rcpp::stop("TransformedProposalKernel::specific_subsample_evaluate_kernel - option for distribution_on_transformed_space not written yet.");
  
  if (this->distribution_on_transformed_space==true)
  {
    //proposed_particle.parameters.deep_overwrite_with_variables_in_argument(this->transform.inverse_transform(proposed_particle.parameters));
    //return this->proposal->subsample_evaluate_kernel(proposed_particle,old_particle) + this->transform.log_abs_inverse_jacobian_determinant(proposed_particle.parameters);
  }
  else
  {
    return this->proposal->subsample_evaluate_kernel(proposed_particle,old_particle);
  }
}

/*
double TransformedProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                     Particle &old_particle,
                                                                     const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
  {
    proposed_particle.parameters.add_parameters(this->transform.inverse_transform(proposed_particle.parameters));
    return this->proposal->subsample_evaluate_kernel(proposed_particle,
                                                     old_particle,
                                                     conditioned_on_parameters) + this->transform.log_abs_inverse_jacobian_determinant(proposed_particle.parameters);
  }
  else
  {
    return this->proposal->subsample_evaluate_kernel(proposed_particle,old_particle,conditioned_on_parameters);
  }
}
*/

Parameters TransformedProposalKernel::simulate(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->simulate(rng,particle));
  else
  {
    Parameters output(this->proposal->simulate(rng,particle));
    output.deep_overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}

/*
Parameters TransformedProposalKernel::simulate(RandomNumberGenerator &rng,
                                               Particle &particle,
                                               const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->simulate(rng,
                                                              particle,
                                                              conditioned_on_parameters));
  else
  {
    Parameters output(this->proposal->simulate(rng,particle,conditioned_on_parameters));
    output.overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}
*/

Parameters TransformedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const Particle &particle) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->subsample_simulate(rng,particle));
  else
  {
    Parameters output(this->proposal->subsample_simulate(rng,particle));
    output.deep_overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}

/*
Parameters TransformedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         Particle &particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->subsample_simulate(rng,particle,
                                                                        conditioned_on_parameters));
  else
  {
    Parameters output(this->proposal->subsample_simulate(rng,particle,conditioned_on_parameters));
    output.overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}
*/

Parameters TransformedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const std::string &variable,
                                                         const Particle &particle) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->subsample_simulate(rng,variable,particle));
  else
  {
    Parameters output(this->proposal->subsample_simulate(rng,variable,particle));
    output.deep_overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}

/*
Parameters TransformedProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const std::string &variable,
                                                         Particle &particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  if (this->distribution_on_transformed_space==true)
    return this->transform.transform(this->proposal->subsample_simulate(rng,
                                                                        variable,
                                                                        particle,
                                                                        conditioned_on_parameters));
  else
  {
    Parameters output(this->proposal->subsample_simulate(rng,variable,particle,conditioned_on_parameters));
    output.overwrite_with_variables_in_argument(this->transform.transform(output));
    return output;
  }
}
*/

arma::mat TransformedProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                              const Particle &proposed_particle,
                                                              const Particle &old_particle)
{
  Rcpp::stop("TransformedProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat TransformedProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                              Particle &proposed_particle,
                                                              Particle &old_particle,
                                                              const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("TransformedProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

//virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
//                                                     Particle &proposed_particle,
//                                                     Particle &old_particle)=0;

arma::mat TransformedProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                        const Particle &proposed_particle,
                                                                        const Particle &old_particle)
{
  Rcpp::stop("TransformedProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat TransformedProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                        Particle &proposed_particle,
                                                                        Particle &old_particle,
                                                                        const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("TransformedProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

void TransformedProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  this->proposal->set_proposal_parameters(proposal_parameters_in);
}

GradientEstimatorOutput* TransformedProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> TransformedProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals = this->proposal->get_proposals();
  proposals.push_back(this);
  return proposals;
}

void TransformedProposalKernel::set_index(Index* index_in)
{
  this->proposal->set_index(index_in);
}

void TransformedProposalKernel::set_index_if_null(Index* index_in)
{
  this->proposal->set_index_if_null(index_in);
}

bool TransformedProposalKernel::can_be_evaluated() const
{
  return proposal->can_be_evaluated();
}
