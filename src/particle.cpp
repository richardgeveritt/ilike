#include "particle.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"

Particle::Particle()
{
}

Particle::Particle(const Parameters &parameters_in)
{
  this->parameters = parameters_in;
}

Particle::~Particle()
{

}

void Particle::simulate()
{

}
