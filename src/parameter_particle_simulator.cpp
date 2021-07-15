#include <iostream>
#include "parameter_particle_simulator.h"
#include "likelihood_estimator.h"
#include "parameters.h"

//Default constructor.
ParameterParticleSimulator::ParameterParticleSimulator()
  :ParticleSimulator()
{
}

ParameterParticleSimulator::ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in,
                                                       const std::vector<LikelihoodEstimator*> &likelihood_estimators_in)
  :ParticleSimulator()
{
  this->simulate_parameters = simulate_parameters_in;
  this->likelihood_estimators = likelihood_estimators_in;
}

//Copy constructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::ParameterParticleSimulator(const ParameterParticleSimulator &another)
  :ParticleSimulator(another)
{
  this->make_copy(another);
}

//Destructor for the ParameterParticleSimulator class.
ParameterParticleSimulator::~ParameterParticleSimulator(void)
{
}

void ParameterParticleSimulator::operator=(const ParameterParticleSimulator &another)
{
  if(this == &another){ //if a==a
    return;
  }

  ParticleSimulator::operator=(another);
  this->make_copy(another);
}

ParticleSimulator* ParameterParticleSimulator::duplicate(void)const
{
  return( new ParameterParticleSimulator(*this));
}

void ParameterParticleSimulator::make_copy(const ParameterParticleSimulator &another)
{
  this->simulate_parameters = another.simulate_parameters;
  this->likelihood_estimators = another.likelihood_estimators;
}

Particle ParameterParticleSimulator::operator()(RandomNumberGenerator &rng)
{
  Parameters simulated_parameters = this->simulate_parameters(rng);
  // Need to have made output before now (means we need to have already made the particle, since this stores the llhd estimator outputs that know how to do the simulation). Then call simulate_auxiliary_variables on this (don't need estimator here in this case).
  
  std::vector<LikelihoodEstimatorOutput*> outputs;
  outputs.reserve(this->likelihood_estimators.size());
  for (std::vector<LikelihoodEstimator*>::iterator i = this->likelihood_estimators.begin();
       i != this->likelihood_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initial_simulate(simulated_parameters));
  }
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, outputs);
}
