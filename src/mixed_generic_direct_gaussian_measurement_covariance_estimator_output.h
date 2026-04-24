#ifndef MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

namespace ilike
{
  /**
   * @file mixed_generic_direct_gaussian_measurement_covariance_estimator_output.h
   * @brief Defines the MixedGenericDirectGaussianMeasurementCovarianceEstimator class.
   *
   * Estimates the mixed generic direct gaussian measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class MixedGenericDirectGaussianMeasurementCovarianceEstimator
   * @brief The mixed generic direct gaussian measurement covariance estimator class.
   */


class MixedGenericDirectGaussianMeasurementCovarianceEstimator;

class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the mixedgenericdirectgaussianmeasurementcovarianceestimatoroutput operation.
   */
  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput();
  /**
   * @brief Performs the ~mixedgenericdirectgaussianmeasurementcovarianceestimatoroutput operation.
   */
  virtual ~MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the mixedgenericdirectgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param estimator_in The estimator.
   */
  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(MixedGenericDirectGaussianMeasurementCovarianceEstimator* estimator_in);
  
  /**
   * @brief Performs the mixedgenericdirectgaussianmeasurementcovarianceestimatoroutput operation.
   *
   * @param another The MixedGenericDirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for MixedGenericDirectGaussianMeasurementCovarianceEstimator.
   *
   * @param another The MixedGenericDirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this MixedGenericDirectGaussianMeasurementCovarianceEstimator object.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  
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
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const arma::mat &inv_sigma_precomp,
                                            double log_det_precomp) const;
  /*
   double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const arma::mat &inv_sigma_precomp,
                                                      double log_det_precomp) const;
  
  /*
   double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
   */
  
  double evaluate_likelihood();
  //double evaluate_likelihood(const Parameters &conditioned_on_parameters);
  /**
   * @brief Performs the subsample evaluate likelihood operation.
   *
   * @return The result.
   */
  double subsample_evaluate_likelihood();
  //double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Returns the estimator.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimator* get_estimator();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  /** @brief The estimator. */
  MixedGenericDirectGaussianMeasurementCovarianceEstimator* estimator;
  
  //Parameters simulated_measurement;
  /** @brief The likelihood measurement state. */
  arma::colvec likelihood_measurement_state;
  /** @brief The prior measurement state. */
  arma::colvec prior_measurement_state;
  /** @brief The likelihood random shift. */
  arma::colvec likelihood_random_shift;
  /** @brief The prior random shift. */
  arma::colvec prior_random_shift;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);
  
  /**
   * @brief Copies the state of another MixedGenericDirectGaussianMeasurementCovarianceEstimator into this object.
   *
   * @param another The MixedGenericDirectGaussianMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
