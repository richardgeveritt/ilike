#include "factor_variables.h"

namespace ilike
{
FactorVariables::FactorVariables()
{
  this->particle = NULL;
}

FactorVariables::FactorVariables(Particle* particle_in)
{
  this->particle = particle_in;
}

FactorVariables::~FactorVariables()
{

}

FactorVariables::FactorVariables(const FactorVariables &another)
{
  this->make_copy(another);
}

void FactorVariables::operator=(const FactorVariables &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void FactorVariables::make_copy(const FactorVariables &another)
{
  this->particle = another.particle;
}

void FactorVariables::set_particle(Particle* particle_in)
{
  this->particle = particle_in;
}

Particle* FactorVariables::get_particle()
{
  return this->particle;
}
}
