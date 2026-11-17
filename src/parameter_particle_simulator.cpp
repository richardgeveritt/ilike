#include <iostream>
#include "parameter_particle_simulator.h"
#include "likelihood_estimator.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"
#include "distributions.h"
#include "independent_proposal_kernel.h"
#include "factors.h"
#include "vector_factor_variables.h"
#include "transform.h"

namespace ilike
{
//Default constructor.
ParameterParticleSimulator::ParameterParticleSimulator()
:ParticleSimulator()
{
  this->proposal.reset();
}

ParameterParticleSimulator::ParameterParticleSimulator(const IndependentProposalKernel* proposal_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
:ParticleSimulator()
{
  this->proposal.reset(proposal_in);
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
  this->proposal.reset();
}

void ParameterParticleSimulator::operator=(const ParameterParticleSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal.reset();
  
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
    this->proposal.reset(another.proposal->independent_proposal_kernel_duplicate());
  else
    this->proposal.reset();
  this->likelihood_estimators = another.likelihood_estimators;
}

Particle ParameterParticleSimulator::simulate(RandomNumberGenerator &rng,
                                              const Factors* factors,
                                              const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for) const
{
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).

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

  // Outputs are created here, with memory managed by Particle hereafter.

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
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).

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
  
  //Parameters simulated_parameters = this->proposal->subsample_independent_simulate(rng,
  //                                                                                 conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).

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
  
  //Parameters simulated_parameters = this->proposal->subsample_independent_simulate(rng,
  //                                                                                 conditioned_on_parameters);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).

  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(this->proposal->subsample_independent_simulate(rng,
                                                                 conditioned_on_parameters),
                  factors,
                  proposals_to_transform_for,
                  proposals_to_find_gradient_for,
                  conditioned_on_parameters,
                  sequencer_parameters);
  
}

void ParameterParticleSimulator::simulate_into(RandomNumberGenerator &rng,
                                               Particle &target,
                                               const Factors* factors,
                                               const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                               const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                               const Parameters &sequencer_parameters) const
{
  // as simulate(...): Particle(simulated, factors, ..., sequencer_parameters)
  target.refill(this->proposal->independent_simulate(rng),
                factors,
                proposals_to_transform_for,
                proposals_to_find_gradient_for,
                sequencer_parameters);
}

void ParameterParticleSimulator::simulate_into(RandomNumberGenerator &rng,
                                               Particle &target,
                                               const Factors* factors,
                                               const std::vector<const ProposalKernel*>* proposals_to_transform_for,
                                               const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for,
                                               const Parameters &conditioned_on_parameters,
                                               const Parameters &sequencer_parameters) const
{
  // as simulate(..., conditioned_on_parameters, sequencer_parameters)
  target.refill(this->proposal->independent_simulate(rng,
                                                     conditioned_on_parameters),
                factors,
                proposals_to_transform_for,
                proposals_to_find_gradient_for,
                conditioned_on_parameters,
                &sequencer_parameters);
}

double ParameterParticleSimulator::evaluate(const Particle &input) const
{
  // could generalise to evaluating the proposals for the likelihoodestimators
  return this->proposal->evaluate_independent_kernel(input.parameters);
}

double ParameterParticleSimulator::subsample_evaluate(const Particle &input) const
{
  return this->proposal->subsample_evaluate_independent_kernel(input.parameters);
}

}
