#ifndef SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H
#define SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "matrix_parameter_estimator.h"

namespace ilike
{
class SampleCovarianceMatrixParameterEstimator : public MatrixParameterEstimator
{
  
public:
  
  SampleCovarianceMatrixParameterEstimator();
  
  virtual ~SampleCovarianceMatrixParameterEstimator();
  
  SampleCovarianceMatrixParameterEstimator(const SampleCovarianceMatrixParameterEstimator &another);
  
  void operator=(const SampleCovarianceMatrixParameterEstimator &another);
  MatrixParameterEstimator* duplicate() const;
  
  void fit(const arma::mat &matrix_points,
           const arma::colvec &normalised_weights);
  
protected:
  
  void make_copy(const SampleCovarianceMatrixParameterEstimator &another);
  
};
}

#endif
