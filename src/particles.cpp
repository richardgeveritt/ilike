#include "particles.h"

//#include <array>

Particles::Particles()
{
}

Particles::Particles(const std::vector<Particle> &particles_in)
{
  this->particles = particles_in;
}

Particles::Particles(const Particles &another)
{
  this->make_copy(another);
}

Particles::~Particles(void)
{
}

void Particles::operator=(const Particles &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Particles::make_copy(const Particles &another)
{
  this->particles = another.particles;
}
