#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"

namespace ilike
{
EnsembleFactorVariables::EnsembleFactorVariables()
{
  this->particle = NULL;
  //this->ensemble_factors = NULL;
}

/*
 EnsembleFactorVariables::EnsembleFactorVariables(EnsembleFactors* ensemble_factors_in)
 {
 this->ensemble_factors = ensemble_factors_in;
 }
 */

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
  //this->ensemble_factors = another.ensemble_factors;
}

/*
 EnsembleFactors* EnsembleFactorVariables::get_ensemble_factors()
 {
 return this->ensemble_factors;
 }
 */

void EnsembleFactorVariables::set_particle(Particle* particle_in)
{
  this->particle = particle_in;
}

Particle* EnsembleFactorVariables::get_particle()
{
  return this->particle;
}
}
