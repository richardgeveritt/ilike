#ifndef ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define ENKGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"

namespace ilike
{
  /**
   * @file enk_gaussian_measurement_covariance_estimator.h
   * @brief Defines the EnKGaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the en k gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class EnKGaussianMeasurementCovarianceEstimator
   * @brief An en k gaussian measurement covariance estimator derived from GaussianMeasurementCovarianceEstimator.
   */


class EnKGaussianMeasurementCovarianceEstimator : public GaussianMeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Default constructor for EnKGaussianMeasurementCovarianceEstimator.
   */
  EnKGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Destructor for EnKGaussianMeasurementCovarianceEstimator.
   */
  virtual ~EnKGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Copy constructor for EnKGaussianMeasurementCovarianceEstimator.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  EnKGaussianMeasurementCovarianceEstimator(const EnKGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for EnKGaussianMeasurementCovarianceEstimator.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const EnKGaussianMeasurementCovarianceEstimator &another);
  /**
   * @brief Creates a deep copy of this EnKGaussianMeasurementCovarianceEstimator object.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  /**
   * @brief Returns the cygivenx.
   *
   * @return The result.
   */
  arma::mat get_Cygivenx() const;
  
  /**
   * @brief Returns the measurement covariance.
   *
   * @return The result.
   */
  arma::mat get_measurement_covariance() const;
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Dhathalf,
                           const arma::mat &P,
                           const arma::mat &Vtranspose,
                           const arma::mat &Yhat,
                           double inverse_incremental_temperature) const;
  
  arma::mat get_sqrt_adjustment(const arma::mat &Sigma,
                                const arma::mat &HSigmaHt,
                                double inverse_incremental_temperature) const;
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  
protected:
  
  /**
   * @brief Performs the initialise measurement covariance estimator operation.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  /**
   * @brief Performs the initialise measurement covariance estimator operation.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the setup measurement variables operation.
   */
  void setup_measurement_variables();
  /**
   * @brief Performs the setup measurement variables operation.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Copies the state of another EnKGaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The EnKGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const EnKGaussianMeasurementCovarianceEstimator &another);
  
};
}

#endif
