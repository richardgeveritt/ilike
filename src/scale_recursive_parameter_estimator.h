#ifndef SCALERECURSIVEPARAMETERESTIMATOR_H
#define SCALERECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "scale.h"
#include "recursive_parameter_estimator.h"

class ScaleRecursiveParameterEstimator : public RecursiveParameterEstimator
{

public:

  ScaleRecursiveParameterEstimator();

  virtual ~ScaleRecursiveParameterEstimator();

  ScaleRecursiveParameterEstimator(const ScaleRecursiveParameterEstimator &another);

  void operator=(const ScaleRecursiveParameterEstimator &another);
  virtual ScaleRecursiveParameterEstimator* scale_duplicate() const=0;
  
  Scale estimated;
  
protected:

  void make_copy(const ScaleRecursiveParameterEstimator &another);

};

#endif
