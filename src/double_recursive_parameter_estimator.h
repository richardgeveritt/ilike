#ifndef MATRIXPARAMETERESTIMATOR_H
#define MATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "recursive_parameter_estimator.h"

namespace ilike
{
class DoubleRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  DoubleRecursiveParameterEstimator();
  DoubleRecursiveParameterEstimator(double initial_value);
  
  virtual ~DoubleRecursiveParameterEstimator();
  
  DoubleRecursiveParameterEstimator(const DoubleRecursiveParameterEstimator &another);
  
  void operator=(const DoubleRecursiveParameterEstimator &another);
  virtual DoubleRecursiveParameterEstimator* double_duplicate() const=0;
  
  double estimated;
  
protected:
  
  void make_copy(const DoubleRecursiveParameterEstimator &another);
  
};
}

#endif
