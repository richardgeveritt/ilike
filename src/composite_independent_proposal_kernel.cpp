#include <iterator>
#include "composite_independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"
#include "particle.h"

CompositeIndependentProposalKernel::CompositeIndependentProposalKernel()
  :IndependentProposalKernel()
{
}

CompositeIndependentProposalKernel::~CompositeIndependentProposalKernel()
{
  for (std::vector<IndependentProposalKernel*>::iterator i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

CompositeIndependentProposalKernel::CompositeIndependentProposalKernel(const std::vector<IndependentProposalKernel*> &all_kernels_in)
:IndependentProposalKernel()
{
  this->all_kernels = all_kernels_in;
}

CompositeIndependentProposalKernel::CompositeIndependentProposalKernel(const CompositeIndependentProposalKernel &another)
  :IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void CompositeIndependentProposalKernel::operator=(const CompositeIndependentProposalKernel &another)
{
  if(this == &another)
    return;
  
  for (std::vector<IndependentProposalKernel*>::iterator i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->all_kernels.clear();

  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* CompositeIndependentProposalKernel::duplicate() const
{
  return( new CompositeIndependentProposalKernel(*this));
}

ProposalKernel* CompositeIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return( new CompositeIndependentProposalKernel(*this));
}

IndependentProposalKernel* CompositeIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return( new CompositeIndependentProposalKernel(*this));
}

void CompositeIndependentProposalKernel::make_copy(const CompositeIndependentProposalKernel &another)
{
  for (std::vector<IndependentProposalKernel*>::const_iterator i=another.all_kernels.begin();
       i!=another.all_kernels.end();
       ++i)
  {
    if (*i!=NULL)
      this->all_kernels.push_back((*i)->independent_proposal_kernel_duplicate());
    else
      this->all_kernels.push_back(NULL);
  }
}

double CompositeIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  double output = 0.0;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    output = output + (*i)->evaluate_independent_kernel(proposed_particle);
  }
  return output;
}

double CompositeIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->evaluate_independent_kernel(proposed_particle);
}

Parameters CompositeIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  Particle output;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->independent_simulate(rng,output.parameters);
  }
  return *output.move_parameters;
}

Parameters CompositeIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                    const Parameters &conditioned_on_parameters) const
{
  Particle output;
  output.parameters.merge(conditioned_on_parameters);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->independent_simulate(rng,output.parameters);
  }
  return *output.move_parameters;
}

Parameters CompositeIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  Particle output;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_independent_simulate(rng,output.parameters);
  }
  return *output.move_parameters;
}

Parameters CompositeIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const Parameters &conditioned_on_parameters) const
{
  Particle output;
  output.parameters.merge(conditioned_on_parameters);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_independent_simulate(rng,output.parameters);
  }
  return *output.move_parameters;
}

Parameters CompositeIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const std::string &variable) const
{
  Particle output;
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_independent_simulate(rng,variable,output.parameters);
  }
  return *output.move_parameters;
}

Parameters CompositeIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                             const std::string &variable,
                                                                             const Parameters &conditioned_on_parameters) const
{
  Particle output;
  output.parameters.merge(conditioned_on_parameters);
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    *output.move_parameters = (*i)->subsample_independent_simulate(rng,variable,output.parameters);
  }
  return *output.move_parameters;
}

arma::mat CompositeIndependentProposalKernel::independent_gradient_of_log(const std::string &variable,
                                                                         const Parameters &proposed_particle)
{
  Rcpp::stop("CompositeIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

arma::mat CompositeIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable,
                                                                                   const Parameters &proposed_particle)
{
  Rcpp::stop("CompositeIndependentProposalKernel::independent_gradient_of_log - not written yet.");
}

void CompositeIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  for (auto i=this->all_kernels.begin();
       i!=this->all_kernels.end();
       ++i)
  {
    (*i)->set_proposal_parameters(proposal_parameters_in);
  }
}
