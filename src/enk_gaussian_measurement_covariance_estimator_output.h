#ifndef ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator_output.h"

namespace ilike
{
  /**
   * @file enk_gaussian_measurement_covariance_estimator_output.h
   * @brief Defines the EnKGaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the en k gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class EnKGaussianMeasurementCovarianceEstimator
   * @brief The en k gaussian measurement covariance estimator class.
   */


class EnKGaussianMeasurementCovarianceEstimator;

class EnKGaussianMeasurementCovarianceEstimatorOutput : public GaussianMeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the enkgaussianmeasurementcovarianceestimatoroutput operation.
   */
  EnKGaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the enkgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param enk_estimator_in The enk estimator.
   */
  EnKGaussianMeasurementCovarianceEstimatorOutput(EnKGaussianMeasurementCovarianceEstimator* enk_estimator_in);
  
  /**
   * @brief Performs the ~enkgaussianmeasurementcovarianceestimatoroutput operation.
   */
  virtual ~EnKGaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the enkgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  EnKGaussianMeasurementCovarianceEstimatorOutput(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for EnKGaussianMeasurementCovarianceEstimator.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this EnKGaussianMeasurementCovarianceEstimator object.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const;
  
  //arma::rowvec get_measurement_state_for_covariance() const;
  //arma::rowvec get_measurement_random_shift();
  
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param parameters The parameters.
   */
  void specific_simulate(const Parameters &parameters);
  /**
   * @brief Performs the subsample specific simulate operation.
   *
   * @param parameters The parameters.
   */
  void subsample_specific_simulate(const Parameters &parameters);
  
  /**
   * @brief Returns the estimator.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimator* get_estimator();
  /**
   * @brief Returns the gaussian estimator.
   *
   * @return The result.
   */
  GaussianMeasurementCovarianceEstimator* get_gaussian_estimator();
  /**
   * @brief Returns the gaussian estimator.
   *
   * @return The result.
   */
  const GaussianMeasurementCovarianceEstimator* get_gaussian_estimator() const;
  
protected:
  
  // not stored here
  /** @brief The enk estimator. */
  EnKGaussianMeasurementCovarianceEstimator* enk_estimator;
  
  /**
   * @brief Copies the state of another EnKGaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const EnKGaussianMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
