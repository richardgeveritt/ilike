#ifndef MATRIXRECURSIVEPARAMETERESTIMATOR_H
#define MATRIXRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "recursive_parameter_estimator.h"

namespace ilike
{
class MatrixRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  MatrixRecursiveParameterEstimator();
  
  virtual ~MatrixRecursiveParameterEstimator();
  
  MatrixRecursiveParameterEstimator(const MatrixRecursiveParameterEstimator &another);
  
  void operator=(const MatrixRecursiveParameterEstimator &another);
  virtual MatrixRecursiveParameterEstimator* matrix_duplicate() const=0;
  
  arma::mat estimated;
  
protected:
  
  void make_copy(const MatrixRecursiveParameterEstimator &another);
  
};
}

#endif
