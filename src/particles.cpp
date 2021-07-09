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

std::ostream& operator<<(std::ostream& os, const Particles &p)
{
  std::vector<Particle>::const_iterator it;

  for (it=p.particles.begin();it!=p.particles.end();++it)
  {
    if (it==p.particles.begin())
      os << *it;
    else
      os << std::endl << *it;
  }

  return os;
}
