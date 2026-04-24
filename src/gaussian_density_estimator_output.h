#ifndef GAUSSIANDENSITYESTIMATOROUTPUT_H
#define GAUSSIANDENSITYESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "density_estimator_output.h"

namespace ilike
{
  /**
   * @file gaussian_density_estimator_output.h
   * @brief Defines the VectorParameterEstimator class.
   *
   * Estimates vector parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class VectorParameterEstimator
   * @brief The vector parameter estimator class.
   */


class VectorParameterEstimator;
class MatrixParameterEstimator;
class GaussianDensityEstimator;

class GaussianDensityEstimatorOutput : public DensityEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the gaussiandensityestimatoroutput operation.
   */
  GaussianDensityEstimatorOutput();
  
  /**
   * @brief Performs the gaussiandensityestimatoroutput operation.
   *
   * @param estimator_in The estimator.
   */
  GaussianDensityEstimatorOutput(GaussianDensityEstimator* estimator_in);
  
  /**
   * @brief Performs the ~gaussiandensityestimatoroutput operation.
   */
  virtual ~GaussianDensityEstimatorOutput();
  
  /**
   * @brief Performs the gaussiandensityestimatoroutput operation.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  GaussianDensityEstimatorOutput(const GaussianDensityEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for VectorParameterEstimator.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void operator=(const GaussianDensityEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this VectorParameterEstimator object.
   *
   * @return The result.
   */
  DensityEstimatorOutput* duplicate() const;
  
  void fit(const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  /**
   * @brief Evaluates.
   *
   * @param point The point.
   *
   * @return The result.
   */
  double evaluate(const Data &point) const;
  
protected:
  
  // not stored heres
  /** @brief The estimator. */
  GaussianDensityEstimator* estimator;
  
  // Stored here.
  /** @brief The mean estimator. */
  VectorParameterEstimator* mean_estimator;
  /** @brief The covariance estimator. */
  MatrixParameterEstimator* covariance_estimator;
  
  /**
   * @brief Copies the state of another VectorParameterEstimator into this object.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void make_copy(const GaussianDensityEstimatorOutput &another);
  
};
}

#endif
