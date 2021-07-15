#ifndef PARTICLES_H
#define PARTICLES_H

//#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
//using namespace RcppParallel;

#include "particle.h"

class Particles
{
public:

  Particles(void);
  Particles(const std::vector<Particle> &particles_in);

  virtual ~Particles(void);

  Particles(const Particles &another);
  void operator=(const Particles &another);

  //virtual List& operator[](int index)=0;
  //virtual List operator[](int index) const=0;

  friend std::ostream& operator<<(std::ostream& os, const Particles &p);

protected:

  void make_copy(const Particles &another);

  std::vector<Particle> particles;

  // Also store anything that needs to be stored in a different structure across all particles.


  // Will loop through these functions to simulate.
  //std::vector<simulatorfiunction> parameters;

  // Stored here. // moved to model_and_algorithm
  //std::vector<LikelihoodEstimator*> likelihood_estimators;

};

#endif
