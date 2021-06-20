//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "parameters.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"

#ifndef PARTICLE_H
#define PARTICLE_H

class Particle
{

protected:

  Parameters parameters;

  // Stored here.
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimators;

  // Pointer, not stored here.
  LikelihoodEstimator* estimator;

public:

  Particle();

  virtual ~Particle();

  // Simulate parameters, then all llhd estimators
  void simulate();

};

#endif
