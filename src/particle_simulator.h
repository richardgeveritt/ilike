#ifndef PARTICLESIMULATOR_H
#define PARTICLESIMULATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "particle.h"

class ParticleSimulator
{
public:

  ParticleSimulator(void);
  ParticleSimulator(const ParticleSimulator &another);
  virtual ~ParticleSimulator(void);

  void operator=(const ParticleSimulator &another);
  virtual ParticleSimulator* duplicate() const=0;

  virtual Particle operator()(RandomNumberGenerator &rng) const=0;

protected:

  void make_copy(const ParticleSimulator &another);

};

#endif
