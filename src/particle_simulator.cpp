#include "particle_simulator.h"

namespace ilike
{
ParticleSimulator::ParticleSimulator()
{
}

ParticleSimulator::ParticleSimulator(const std::string &resample_variable_name_in)
{
  this->resample_variable_name = resample_variable_name_in;
}

ParticleSimulator::ParticleSimulator(const ParticleSimulator &another)
{
  this->make_copy(another);
}

ParticleSimulator::~ParticleSimulator()
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
  this->resample_variable_name = another.resample_variable_name;
}
}
