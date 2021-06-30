#include "parameter_particle_simulator.h"

//Default constructor.
ParameterParticleSimulator::ParameterParticleSimulator(SimulateDistributionPtr simulate_parameters_in)
  :ParticleSimulator()
{
  this->simulate_parameters = simulate_parameters_in;
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
}

Particle ParameterParticleSimulator::operator()(RandomNumberGenerator &rng) const
{
  Parameters simulated_parameters = this->simulate_parameters(rng);
  return Particle(simulated_parameters);
}
