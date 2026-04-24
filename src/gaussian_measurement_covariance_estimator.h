#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"

namespace ilike
{
  /**
   * @file gaussian_measurement_covariance_estimator.h
   * @brief Defines the GaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class GaussianMeasurementCovarianceEstimator
   * @brief A gaussian measurement covariance estimator derived from MeasurementCovarianceEstimator.
   */


class GaussianMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Default constructor for GaussianMeasurementCovarianceEstimator.
   */
  GaussianMeasurementCovarianceEstimator();
  
  GaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                         size_t* seed_in,
                                         Data* data_in,
                                         std::shared_ptr<Transform> inverse_transform_in,
                                         std::shared_ptr<Transform> summary_statistics_in,
                                         const std::vector<std::string> &measurement_variables_in);
  
  /**
   * @brief Destructor for GaussianMeasurementCovarianceEstimator.
   */
  virtual ~GaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Copy constructor for GaussianMeasurementCovarianceEstimator.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  GaussianMeasurementCovarianceEstimator(const GaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for GaussianMeasurementCovarianceEstimator.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const GaussianMeasurementCovarianceEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  virtual GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const=0;
  
  /**
   * @brief Performs the need cxx operation.
   *
   * @return The result.
   */
  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature) const;
  
  /**
   * @brief Returns the measurement covariance.
   *
   * @return The result.
   */
  virtual arma::mat get_measurement_covariance() const=0;
  
  /**
   * @brief Performs the change data operation.
   */
  void change_data();
  /**
   * @brief Performs the change data operation.
   *
   * @param new_data The new data.
   */
  void change_data(std::shared_ptr<Data> new_data);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                      arma::mat &inv_sigma_precomp,
                                      double &log_det_precomp);
  
protected:
  
  /**
   * @brief Copies the state of another GaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const GaussianMeasurementCovarianceEstimator &another);
  
};
}

#endif
