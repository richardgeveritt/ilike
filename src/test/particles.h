#include <RcppParallel.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace RcppParallel;

#include "parameters.h"
#include "likelihood_estimator.h"

#ifndef PARTICLES_H
#define PARTICLES_H

class Particles
{
public:

  Particles(void);
  Particles(const Particles &another);
  virtual ~Particles(void);

  void operator=(const Particles &another);

  //virtual List& operator[](int index)=0;
  //virtual List operator[](int index) const=0;

protected:

  void make_copy(const Particles &another);

  std::vector<Parameters> particles;

  // Also store anything that needs to be stored in a different structure across all particles.


  // Will loop through these functions to simulate.
  //std::vector<simulatorfiunction> parameters;

  // Stored here.
  std::vector<LikelihoodEstimator*> likelihood_estimators;

};

#endif
