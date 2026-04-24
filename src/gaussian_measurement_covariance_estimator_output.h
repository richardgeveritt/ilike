#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

namespace ilike
{
  /**
   * @file gaussian_measurement_covariance_estimator_output.h
   * @brief Defines the GaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class GaussianMeasurementCovarianceEstimator
   * @brief The gaussian measurement covariance estimator class.
   */


class GaussianMeasurementCovarianceEstimator;

class GaussianMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the gaussianmeasurementcovarianceestimatoroutput operation.
   */
  GaussianMeasurementCovarianceEstimatorOutput();
  /**
   * @brief Performs the ~gaussianmeasurementcovarianceestimatoroutput operation.
   */
  virtual ~GaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the gaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  GaussianMeasurementCovarianceEstimatorOutput(const GaussianMeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for GaussianMeasurementCovarianceEstimator.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const GaussianMeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  virtual GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const=0;
  
  /**
   * @brief Returns the measurement state for covariance.
   *
   * @return The result.
   */
  arma::rowvec get_measurement_state_for_covariance() const;
  
  /**
   * @brief Returns the shift.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  
  /**
   * @brief Returns the deterministic shift.
   *
   * @return The result.
   */
  arma::colvec get_deterministic_shift() const;
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature);
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const arma::mat &inv_sigma_precomp,
                                            double log_det_precomp) const;
  
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const arma::mat &inv_sigma_precomp,
                                                      double log_det_precomp) const;
  
  /*
   double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  double evaluate_likelihood();
  /*
   double evaluate_likelihood(const Parameters &conditioned_on_parameters);
   */
  double subsample_evaluate_likelihood();
  /*
   double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);
   */
  
  //virtual arma::mat get_measurement_covariance()=0;
  
  /**
   * @brief Returns the gaussian estimator.
   *
   * @return The result.
   */
  virtual GaussianMeasurementCovarianceEstimator* get_gaussian_estimator()=0;
  /**
   * @brief Returns the gaussian estimator.
   *
   * @return The result.
   */
  virtual const GaussianMeasurementCovarianceEstimator* get_gaussian_estimator() const=0;
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  /** @brief The measurement state. */
  arma::colvec measurement_state;
  /** @brief The random shift. */
  arma::colvec random_shift;
  
  //arma::mat measurement_noise;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);
  
  /**
   * @brief Copies the state of another GaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The GaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const GaussianMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
