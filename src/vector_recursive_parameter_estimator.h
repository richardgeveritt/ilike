#ifndef VECTORRECURSIVEPARAMETERESTIMATOR_H
#define VECTORRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include <boost/unordered_map.hpp>
#include "recursive_parameter_estimator.h"

namespace ilike
{
class VectorRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  VectorRecursiveParameterEstimator();
  
  virtual ~VectorRecursiveParameterEstimator();
  
  VectorRecursiveParameterEstimator(const VectorRecursiveParameterEstimator &another);
  
  void operator=(const VectorRecursiveParameterEstimator &another);
  virtual VectorRecursiveParameterEstimator* vector_duplicate() const=0;
  
  arma::colvec estimated;
  
protected:
  
  void make_copy(const VectorRecursiveParameterEstimator &another);
  
};
}

#endif
