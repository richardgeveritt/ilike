#include <iostream>
#include "parameter_particle_simulator.h"
#include "likelihood_estimator.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"
#include "distributions.h"
#include "independent_proposal_kernel.h"
//#include "custom_independent_proposal_kernel.h"
#include "factors.h"
#include "vector_factor_variables.h"
#include "transform.h"

//Default constructor.
ParameterParticleSimulator::ParameterParticleSimulator()
  :ParticleSimulator()
{
  this->proposal = NULL;
}

/*
ParameterParticleSimulator::ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                                                       const std::string &resample_variable_name_in)
  :ParticleSimulator(resample_variable_name_in)
{
  this->simulate_parameters = simulate_parameters_in;
  this->likelihood_estimators = likelihood_estimators_in;
}
*/

ParameterParticleSimulator::ParameterParticleSimulator(const IndependentProposalKernel* proposal_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:ParticleSimulator()
{
  this->proposal = proposal_in;
  this->likelihood_estimators = likelihood_estimators_in;
}

//Copy constructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::ParameterParticleSimulator(const ParameterParticleSimulator &another)
  :ParticleSimulator(another)
{
  this->make_copy(another);
}

//Destructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::~ParameterParticleSimulator()
{
  if (this->proposal!=NULL)
    delete this->proposal;
}

void ParameterParticleSimulator::operator=(const ParameterParticleSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  this->likelihood_estimators.clear();

  ParticleSimulator::operator=(another);
  this->make_copy(another);
}

ParticleSimulator* ParameterParticleSimulator::duplicate() const
{
  return( new ParameterParticleSimulator(*this));
}

void ParameterParticleSimulator::make_copy(const ParameterParticleSimulator &another)
{
  if (another.proposal!=NULL)
    this->proposal = another.proposal->independent_proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  this->likelihood_estimators = another.likelihood_estimators;
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              const Factors* factors,
                                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const
{
  //new_particle->setup(this->proposal->independent_simulate(rng),
  //                    factors);
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->independent_simulate(rng),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for);
  
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              const Factors* factors,
                                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                              const Parameters &sequencer_parameters) const
{
  //new_particle->setup(this->proposal->independent_simulate(rng,
  //                                                         conditioned_on_parameters),
  //                    factors,
  //                    conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters,
  //                                                                                 conditioned_on_parameters);
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng,
  //                                                                       conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  
  //this->proposal->independent_simulate(rng);
  
  //return Particle();
  
  return Particle(this->proposal->independent_simulate(rng),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for,
                  sequencer_parameters);
  
  
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              const Factors* factors,
                                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                              const Parameters &conditioned_on_parameters,
                                              const Parameters &sequencer_parameters) const
{
  //new_particle->setup(this->proposal->independent_simulate(rng,
  //                                                         conditioned_on_parameters),
  //                    factors,
  //                    conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters,
  //                                                                                 conditioned_on_parameters);
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng,
  //                                                                       conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->independent_simulate(rng,
                                                       conditioned_on_parameters),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for,
                  conditioned_on_parameters,
                  sequencer_parameters);
  
}

Particle ParameterParticleSimulator::subsample_simulate(RandomNumberGenerator &rng,
                                                        const Factors* factors,
                                                        const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const
{
  //new_particle->setup(this->proposal->independent_simulate(rng),
  //                    factors);
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->simulate_factor_variables(simulated_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->subsample_independent_simulate(rng),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for);
  
}

Particle ParameterParticleSimulator::subsample_simulate(RandomNumberGenerator &rng,
                                                        const Factors* factors,
                                                        const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                                        const Parameters &sequencer_parameters) const
{
  //new_particle->setup(this->proposal->subsample_independent_simulate(rng,
  //                                                                   conditioned_on_parameters),
  //                    factors,
  //                    conditioned_on_parameters);
  
  //Parameters simulated_parameters = this->proposal->subsample_independent_simulate(rng,
  //                                                                                 conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->subsample_simulate_factor_variables(simulated_parameters,
  //                                                                                           conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->subsample_independent_simulate(rng),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for,
                  sequencer_parameters);
  
}

Particle ParameterParticleSimulator::subsample_simulate(RandomNumberGenerator &rng,
                                                        const Factors* factors,
                                                        const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                                        const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                                        const Parameters &conditioned_on_parameters,
                                                        const Parameters &sequencer_parameters) const
{
  //new_particle->setup(this->proposal->subsample_independent_simulate(rng,
  //                                                                   conditioned_on_parameters),
  //                    factors,
  //                    conditioned_on_parameters);
  
  //Parameters simulated_parameters = this->proposal->subsample_independent_simulate(rng,
  //                                                                                 conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  //if (!conditioned_on_parameters.is_empty())
  //  simulated_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  //FactorVariables* simulated_factor_variables = factors->subsample_simulate_factor_variables(simulated_parameters,
  //                                                                                           conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->subsample_independent_simulate(rng,
                                                                 conditioned_on_parameters),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for,
                  conditioned_on_parameters,
                  sequencer_parameters);
  
}

/*
void ParameterParticleSimulator::simulate_and_transform(RandomNumberGenerator &rng,
                                                        Particle* new_particle,
                                                        Factors* factors,
                                                        Transform* transform,
                                                        bool store_raw) const
{
  new_particle->parameters = this->proposal->independent_simulate(rng);
  if (!store_raw)
  {
    new_particle->parameters = transform->transform(new_particle->parameters);
  }
  else
  {
    new_particle->parameters.add_parameters_overwrite(transform->transform(new_particle->parameters));
  }
  new_particle->setup(factors);
}

void ParameterParticleSimulator::simulate_and_transform(RandomNumberGenerator &rng,
                                                        Particle* new_particle,
                                                        Factors* factors,
                                                        Transform* transform,
                                                        bool store_raw,
                                                        const Parameters &conditioned_on_parameters) const
{
  new_particle->parameters = this->proposal->independent_simulate(rng,
                                                                  conditioned_on_parameters);
  if (!store_raw)
  {
    new_particle->parameters = transform->transform(new_particle->parameters);
  }
  else
  {
    new_particle->parameters.add_parameters_overwrite(transform->transform(new_particle->parameters));
  }
  new_particle->setup(factors,
                      conditioned_on_parameters);
}

void ParameterParticleSimulator::subsample_simulate_and_transform(RandomNumberGenerator &rng,
                                                                  Particle* new_particle,
                                                                  Factors* factors,
                                                                  Transform* transform,
                                                                  bool store_raw,
                                                                  const Parameters &conditioned_on_parameters) const
{
  new_particle->parameters = this->proposal->subsample_independent_simulate(rng,
                                                                            conditioned_on_parameters);
  if (!store_raw)
  {
    new_particle->parameters = transform->transform(new_particle->parameters);
  }
  else
  {
    new_particle->parameters.add_parameters_overwrite(transform->transform(new_particle->parameters));
  }
  new_particle->setup(factors,
                      conditioned_on_parameters);
}
*/

double ParameterParticleSimulator::evaluate(const Particle &input) const
{
  // could generalise to evaluating the proposals for the likelihoodestimators
  return this->proposal->evaluate_independent_kernel(input.parameters);
}

double ParameterParticleSimulator::subsample_evaluate(const Particle &input) const
{
  return this->proposal->subsample_evaluate_independent_kernel(input.parameters);
}

/*
double ParameterParticleSimulator::evaluate(Particle &input,
                                            const Parameters &conditioned_on_parameters) const
{
  // could generalise to evaluating the proposals for the likelihoodestimators
  if (input.parameters.is_empty())
  {
    return this->proposal->evaluate_independent_kernel(conditioned_on_parameters);
  }
  
  if (conditioned_on_parameters.is_empty())
  {
    return this->proposal->evaluate_independent_kernel(input.parameters);
  }
  
  return this->proposal->evaluate_independent_kernel(input.parameters.merge(conditioned_on_parameters));
}

double ParameterParticleSimulator::subsample_evaluate(Particle &input,
                                                      const Parameters &conditioned_on_parameters) const
{
  if (input.parameters.is_empty())
  {
    return this->proposal->subsample_evaluate_independent_kernel(conditioned_on_parameters);
  }
  
  if (conditioned_on_parameters.is_empty())
  {
    return this->proposal->subsample_evaluate_independent_kernel(input.parameters);
  }
  
  return this->proposal->subsample_evaluate_independent_kernel(input.parameters.merge(conditioned_on_parameters));
}
*/
