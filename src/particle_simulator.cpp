#include "particle_simulator.h"

ParticleSimulator::ParticleSimulator(void)
{
}

ParticleSimulator::ParticleSimulator(const ParticleSimulator &another)
{
  this->make_copy(another);
}

ParticleSimulator::~ParticleSimulator(void)
{
}

void ParticleSimulator::operator=(const ParticleSimulator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void ParticleSimulator::make_copy(const ParticleSimulator &another)
{
  //Does nothing since no member variables to copy.
}
