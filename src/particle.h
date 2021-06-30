#ifndef PARTICLE_H
#define PARTICLE_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;
#include <vector>
#include "parameters.h"

class ModelAndAlgorithm;
class LikelihoodEstimatorOutput;

class Particle
{

public:

  Particle();
  Particle(const Parameters &parameters_in);

  virtual ~Particle();

  // Simulate parameters, then all llhd estimators
  void simulate();

protected:

  Parameters parameters;

  // Stored here.
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;

  // Pointer, not stored here.
  const ModelAndAlgorithm* model_and_algorithm;

};

#endif
