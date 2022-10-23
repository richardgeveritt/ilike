#ifndef DOUBLESCALERECURSIVEPARAMETERESTIMATOR_H
#define DOUBLESCALERECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "particle.h"
#include "scale_recursive_parameter_estimator.h"

class DoubleRecursiveParameterEstimator;

class DoubleScaleRecursiveParameterEstimator : public ScaleRecursiveParameterEstimator
{

public:

  DoubleScaleRecursiveParameterEstimator();

  virtual ~DoubleScaleRecursiveParameterEstimator();

  DoubleScaleRecursiveParameterEstimator(const DoubleScaleRecursiveParameterEstimator &another);

  void operator=(const DoubleScaleRecursiveParameterEstimator &another);
  RecursiveParameterEstimator* duplicate() const;
  ScaleRecursiveParameterEstimator* scale_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  // stored here
  DoubleRecursiveParameterEstimator* recursive_estimator;

  void make_copy(const DoubleScaleRecursiveParameterEstimator &another);

};

#endif
