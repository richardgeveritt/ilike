#ifndef PARTICLE_H
#define PARTICLE_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;
#include <vector>
#include <iostream>
#include "parameters.h"

class ModelAndAlgorithm;
class LikelihoodEstimatorOutput;

class Particle
{

public:

  Particle();
  Particle(const Parameters &parameters_in);
  Particle(const Parameters &parameters_in,
           const std::vector<LikelihoodEstimatorOutput*> &outputs_in);
  virtual ~Particle();

  Particle(const Particle &another);
  void operator=(const Particle &another);

  friend std::ostream& operator<<(std::ostream& os, const Particle &p);

protected:

  void make_copy(const Particle &another);

  Parameters parameters;

  // Stored here.
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;

  // Pointer, not stored here.
  //const ModelAndAlgorithm* model_and_algorithm;

};

#endif
