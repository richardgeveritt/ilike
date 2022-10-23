#ifndef SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H
#define SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "matrix_parameter_estimator.h"

class SampleCovarianceMatrixParameterEstimator : public MatrixParameterEstimator
{

public:

  SampleCovarianceMatrixParameterEstimator();

  virtual ~SampleCovarianceMatrixParameterEstimator();

  SampleCovarianceMatrixParameterEstimator(const SampleCovarianceMatrixParameterEstimator &another);

  void operator=(const SampleCovarianceMatrixParameterEstimator &another);
  MatrixParameterEstimator* duplicate() const;

  void fit(const std::vector<Parameters> &points,
           arma::colvec normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<Parameters> &points,
           arma::colvec normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<MoveOutput*> &points,
           arma::colvec normalised_log_weights);
  
protected:

  void make_copy(const SampleCovarianceMatrixParameterEstimator &another);

};

#endif
