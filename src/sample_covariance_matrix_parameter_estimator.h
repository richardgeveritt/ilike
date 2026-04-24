#ifndef SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H
#define SAMPLECOVARIANCEMATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "matrix_parameter_estimator.h"

namespace ilike
{
  /**
   * @file sample_covariance_matrix_parameter_estimator.h
   * @brief Defines the SampleCovarianceMatrixParameterEstimator class.
   *
   * Estimates sample covariance matrix parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class SampleCovarianceMatrixParameterEstimator
   * @brief A sample covariance matrix parameter estimator derived from MatrixParameterEstimator.
   */


class SampleCovarianceMatrixParameterEstimator : public MatrixParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for SampleCovarianceMatrixParameterEstimator.
   */
  SampleCovarianceMatrixParameterEstimator();
  
  /**
   * @brief Destructor for SampleCovarianceMatrixParameterEstimator.
   */
  virtual ~SampleCovarianceMatrixParameterEstimator();
  
  /**
   * @brief Copy constructor for SampleCovarianceMatrixParameterEstimator.
   *
   * @param another The SampleCovarianceMatrixParameterEstimator instance to copy from.
   */
  SampleCovarianceMatrixParameterEstimator(const SampleCovarianceMatrixParameterEstimator &another);
  
  /**
   * @brief Assignment operator for SampleCovarianceMatrixParameterEstimator.
   *
   * @param another The SampleCovarianceMatrixParameterEstimator instance to copy from.
   */
  void operator=(const SampleCovarianceMatrixParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this SampleCovarianceMatrixParameterEstimator object.
   *
   * @return The result.
   */
  MatrixParameterEstimator* duplicate() const;
  
  void fit(const arma::mat &matrix_points,
           const arma::colvec &normalised_weights);
  
protected:
  
  /**
   * @brief Copies the state of another SampleCovarianceMatrixParameterEstimator into this object.
   *
   * @param another The SampleCovarianceMatrixParameterEstimator instance to copy from.
   */
  void make_copy(const SampleCovarianceMatrixParameterEstimator &another);
  
};
}

#endif
