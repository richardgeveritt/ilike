//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "batch_simulator.h"

#ifndef MODELANDALGORITHM_H
#define MODELANDALGORITHM_H

class ModelAndAlgorithm
{

public:

  std::vector<std::string> is_simulate_methods;

  BatchSimulator* simulate_priors;
  //std::vector<> simulate_for_likelihoods;

  ModelAndAlgorithm();
  ModelAndAlgorithm(const ModelAndAlgorithm &another);

  virtual ~ModelAndAlgorithm();

  void operator=(const ModelAndAlgorithm &another);
  ModelAndAlgorithm* duplicate() const;

protected:

  void make_copy(const ModelAndAlgorithm &another);

};

#endif
