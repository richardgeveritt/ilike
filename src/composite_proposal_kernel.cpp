#include "composite_proposal_kernel.h"
#include "mcmc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "smc_adaptor.h"
#include "transform.h"

CompositeProposalKernel::CompositeProposalKernel()
  :ProposalKernel()
{
}

CompositeProposalKernel::~CompositeProposalKernel()
{
  for (std::vector<ProposalKernel*>::iterator i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

CompositeProposalKernel::CompositeProposalKernel(const std::vector<ProposalKernel*> &all_kernels_in)
:ProposalKernel()
{
  this->all_kernels = all_kernels_in;
}

CompositeProposalKernel::CompositeProposalKernel(const CompositeProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void CompositeProposalKernel::operator=(const CompositeProposalKernel &another)
{
  if(this == &another)
    return;
  
  for (std::vector<ProposalKernel*>::iterator i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->all_kernels.clear();

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CompositeProposalKernel::duplicate() const
{
  return( new CompositeProposalKernel(*this));
}

ProposalKernel* CompositeProposalKernel::proposal_kernel_duplicate() const
{
  return( new CompositeProposalKernel(*this));
}

void CompositeProposalKernel::make_copy(const CompositeProposalKernel &another)
{
  for (std::vector<ProposalKernel*>::const_iterator i=another.all_kernels.begin();
       i!=another.all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      this->all_kernels.push_back((*i)->proposal_kernel_duplicate());
    else
      this->all_kernels.push_back(NULL);
  }
}

double CompositeProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                         const Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    output = output + (*i)->evaluate_kernel(proposed_particle,
                                            old_particle);
  }
  return output;
}

/*
double CompositeProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                         Particle &old_particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  double output = 0.0;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    output = output + (*i)->evaluate_kernel(proposed_particle,
                                            old_particle,
                                            conditioned_on_parameters);
  }
  return output;
}
*/

double CompositeProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                   const Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    output = output + (*i)->subsample_evaluate_kernel(proposed_particle,
                                                      old_particle);
  }
  return output;
}

/*
double CompositeProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                   Particle &old_particle,
                                                                   const Parameters &conditioned_on_parameters) const
{
  double output = 0.0;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    output = output + (*i)->subsample_evaluate_kernel(proposed_particle,
                                                      old_particle,
                                                      conditioned_on_parameters);
  }
  return output;
}
*/

Parameters CompositeProposalKernel::simulate(RandomNumberGenerator &rng,
                                             const Particle &particle) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    Parameters proposed_parameters_in_transformed_space = (*i)->simulate(rng,output);
    output.parameters.deep_overwrite_with_variables_in_argument((*i)->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  return output.parameters;
}

/*
Parameters CompositeProposalKernel::simulate(RandomNumberGenerator &rng,
                                             Particle &particle,
                                             const Parameters &conditioned_on_parameters) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->simulate(rng,output,conditioned_on_parameters);
  }
  return *output.move_parameters;
}
*/

Parameters CompositeProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                       const Particle &particle) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    Parameters proposed_parameters_in_transformed_space = (*i)->subsample_simulate(rng,output);
    output.parameters.deep_overwrite_with_variables_in_argument((*i)->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  return output.parameters;
}

/*
Parameters CompositeProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                       Particle &particle,
                                                       const Parameters &conditioned_on_parameters) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_simulate(rng,output,conditioned_on_parameters);
  }
  return *output.move_parameters;
}
*/

Parameters CompositeProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                       const std::string &variable,
                                                       const Particle &particle) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    Parameters proposed_parameters_in_transformed_space = (*i)->subsample_simulate(rng,variable,output);
    output.parameters.deep_overwrite_with_variables_in_argument((*i)->transform->inverse_transform(proposed_parameters_in_transformed_space));
  }
  return output.parameters;
}

/*
Parameters CompositeProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                       const std::string &variable,
                                                       Particle &particle,
                                                       const Parameters &conditioned_on_parameters) const
{
  Particle output(particle);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_simulate(rng,
                                                       variable,
                                                       output,
                                                       conditioned_on_parameters);
  }
  return *output.move_parameters;
}
*/

arma::mat CompositeProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                            const Particle &proposed_particle,
                                                            const Particle &old_particle)
{
  Rcpp::stop("CompositeProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat CompositeProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                      const Particle &proposed_particle,
                                                                      const Particle &old_particle)
{
  Rcpp::stop("CompositeProposalKernel::specific_subsample_gradient_of_log - not written yet.");
}

void CompositeProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}

GradientEstimatorOutput* CompositeProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> CompositeProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> all_proposals;
  
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    std::vector<const ProposalKernel*> next_proposals = (*i)->get_proposals();
    if (all_proposals.size()==0)
    {
      all_proposals.insert(all_proposals.end(), next_proposals.begin(), next_proposals.end());
    }
  }
  
  all_proposals.push_back(this);
  
  return all_proposals;
}

void CompositeProposalKernel::set_index(Index* index_in)
{
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    (*i)->set_index(index_in);
  }
}

void CompositeProposalKernel::set_index_if_null(Index* index_in)
{
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    (*i)->set_index_if_null(index_in);
  }
}

/*
arma::mat CompositeProposalKernel::specific_gradient_of_log(const std::string &variable,
                                           Particle &proposed_particle,
                                           Particle &old_particle,
                                           const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CompositeProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

//virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
//                                                     Particle &proposed_particle,
//                                                     Particle &old_particle)=0;

/*
arma::mat CompositeProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                     Particle &proposed_particle,
                                                     Particle &old_particle,
                                                     const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("CompositeProposalKernel::specific_gradient_of_log - not written yet.");
}
*/
