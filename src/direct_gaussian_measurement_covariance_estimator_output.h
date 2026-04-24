#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

//#include <RcppArmadillo.h>
//using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator_output.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file direct_gaussian_measurement_covariance_estimator_output.h
   * @brief Defines the DirectGaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the direct gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class DirectGaussianMeasurementCovarianceEstimator
   * @brief The direct gaussian measurement covariance estimator class.
   */


class DirectGaussianMeasurementCovarianceEstimator;

class DirectGaussianMeasurementCovarianceEstimatorOutput : public GaussianMeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the directgaussianmeasurementcovarianceestimatoroutput operation.
   */
  DirectGaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the directgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param direct_estimator_in The direct estimator.
   */
  DirectGaussianMeasurementCovarianceEstimatorOutput(DirectGaussianMeasurementCovarianceEstimator* direct_estimator_in);
  
  /**
   * @brief Performs the ~directgaussianmeasurementcovarianceestimatoroutput operation.
   */
  virtual ~DirectGaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the directgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  DirectGaussianMeasurementCovarianceEstimatorOutput(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for DirectGaussianMeasurementCovarianceEstimator.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this DirectGaussianMeasurementCovarianceEstimator object.
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
  
  //arma::rowvec get_measurement_random_shift();
  
  //arma::mat get_measurement_covariance() const;
  
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
  /** @brief The direct estimator. */
  DirectGaussianMeasurementCovarianceEstimator* direct_estimator;
  
  /**
   * @brief Copies the state of another DirectGaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
