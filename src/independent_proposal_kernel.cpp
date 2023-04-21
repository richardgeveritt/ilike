#include "independent_proposal_kernel.h"
#include "smc_adaptor.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "mcmc_adaptor.h"

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
  //this->proposal_evaluate = another.proposal_evaluate;
  //this->proposal_simulate = another.proposal_simulate;
}

double IndependentProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                           Particle &old_particle) const
{
  return this->evaluate_independent_kernel(*proposed_particle.move_parameters);
}

/*
double IndependentProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                           Particle &old_particle,
                                                           const Parameters &conditioned_on_parameters) const
{
  return this->evaluate_independent_kernel(proposed_particle.move_parameters->merge(conditioned_on_parameters));
}
*/

double IndependentProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                     Particle &old_particle) const
{
  return this->subsample_evaluate_independent_kernel(*proposed_particle.move_parameters);
}

/*
double IndependentProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                     Particle &old_particle,
                                                                     const Parameters &conditioned_on_parameters) const
{
  return this->subsample_evaluate_independent_kernel(proposed_particle.move_parameters->merge(conditioned_on_parameters));
}
*/

arma::mat IndependentProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                              Particle &proposed_particle,
                                                              Particle &old_particle)
{
  return this->independent_gradient_of_log(variable,
                                           *proposed_particle.move_parameters);
}

/*
arma::mat IndependentProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                              Particle &proposed_particle,
                                                              Particle &old_particle,
                                                              const Parameters &conditioned_on_parameters)
{
  return this->independent_gradient_of_log(variable,
                                           proposed_particle.move_parameters->merge(conditioned_on_parameters));
}
*/

arma::mat IndependentProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                        Particle &proposed_particle,
                                                                        Particle &old_particle)
{
  return this->subsample_independent_gradient_of_log(variable,
                                                     *proposed_particle.move_parameters);
}

/*
arma::mat IndependentProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                        Particle &proposed_particle,
                                                                        Particle &old_particle,
                                                                        const Parameters &conditioned_on_parameters)
{
  return this->subsample_independent_gradient_of_log(variable,
                                                     proposed_particle.move_parameters->merge(conditioned_on_parameters));
}
*/

Parameters IndependentProposalKernel::simulate(RandomNumberGenerator &rng,
                                               Particle &particle) const
{
  Parameters output = particle.parameters;
  output.deep_overwrite_with_variables_in_argument(this->independent_simulate(rng));
  return output;
}

/*
Parameters IndependentProposalKernel::simulate(RandomNumberGenerator &rng,
                                               Particle &particle,
                                               const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng,
                                    conditioned_on_parameters);
}
*/

Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         Particle &particle) const
{
  Parameters output = particle.parameters;
  output.deep_overwrite_with_variables_in_argument(this->subsample_independent_simulate(rng));
  return output;
}

/*
Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         Particle &particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  return this->subsample_independent_simulate(rng,
                                              conditioned_on_parameters);
}
*/

Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const std::string &variable,
                                                         Particle &particle) const
{
  Parameters output = particle.parameters;
  output.deep_overwrite_with_variables_in_argument(this->subsample_independent_simulate(rng,variable));
  return output;
}

/*
Parameters IndependentProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                         const std::string &variable,
                                                         Particle &particle,
                                                         const Parameters &conditioned_on_parameters) const
{
  return this->subsample_independent_simulate(rng,
                                              variable,
                                              conditioned_on_parameters);
}
*/
