#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"

namespace ilike
{
EnsembleFactorVariables::EnsembleFactorVariables()
{
  this->particle = NULL;
}

EnsembleFactorVariables::~EnsembleFactorVariables()
{
  
}

EnsembleFactorVariables::EnsembleFactorVariables(const EnsembleFactorVariables &another)
{
  this->make_copy(another);
}

void EnsembleFactorVariables::operator=(const EnsembleFactorVariables &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void EnsembleFactorVariables::make_copy(const EnsembleFactorVariables &another)
{
  this->particle = another.particle;
}

void EnsembleFactorVariables::set_particle(Particle* particle_in)
{
  this->particle = particle_in;
}

Particle* EnsembleFactorVariables::get_particle()
{
  return this->particle;
}
}
