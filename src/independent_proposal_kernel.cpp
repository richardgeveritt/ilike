#include "independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

namespace ilike
{
IndependentProposalKernel::IndependentProposalKernel()
:ProposalKernel()
{
}

IndependentProposalKernel::~IndependentProposalKernel()
{
}

IndependentProposalKernel::IndependentProposalKernel(const IndependentProposalKernel &another)
:ProposalKernel(another)
{
  this->make_copy(another);
}

void IndependentProposalKernel::operator=(const IndependentProposalKernel &another)
{
  if(this == &another)
    return;
  
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

void IndependentProposalKernel::make_copy(const IndependentProposalKernel &another)
{
}

double IndependentProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                           const Particle &old_particle) const
{
  return this->evaluate_independent_kernel(proposed_particle.get_transformed_parameters(this));
}

double IndependentProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                     const Particle &old_particle) const
{
  return this->subsample_evaluate_independent_kernel(proposed_particle.get_transformed_parameters(this));
}

arma::mat IndependentProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                              const Particle &proposed_particle,
                                                              const Particle &old_particle)
{
  return this->independent_gradient_of_log(variable,
                                           proposed_particle.get_transformed_parameters(this));
}

arma::mat IndependentProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                        const Particle &proposed_particle,
                                                                        const Particle &old_particle)
{
  return this->subsample_independent_gradient_of_log(variable,
                                                     proposed_particle.get_transformed_parameters(this));
}

// in a move (independent_mh): drawn for the blocks at the current point, and
// each variable given the shape it has there (a stacked name split into its
// variables, C205)
Parameters IndependentProposalKernel::simulate(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  Parameters drawn = this->independent_simulate_like(rng, particle.parameters);
  Parameters output;
  for (auto i = drawn.vector_begin(); i != drawn.vector_end(); ++i)
    output.set_stacked(i->first, arma::vectorise(*i->second.first), particle.parameters);
  for (auto i = drawn.any_begin(); i != drawn.any_end(); ++i)
    output(i->first) = *i->second.first;
  return output;
}

Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const Particle &particle) const
{
  return this->subsample_independent_simulate(rng);
}

Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const std::string &variable,
                                                         const Particle &particle) const
{
  return this->subsample_independent_simulate(rng, variable);
}

}
