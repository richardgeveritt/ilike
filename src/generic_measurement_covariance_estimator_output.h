#ifndef GENERICMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define GENERICMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

namespace ilike
{
  /**
   * @file generic_measurement_covariance_estimator_output.h
   * @brief Defines the GenericMeasurementCovarianceEstimator class.
   *
   * Estimates the generic measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class GenericMeasurementCovarianceEstimator
   * @brief The generic measurement covariance estimator class.
   */


class GenericMeasurementCovarianceEstimator;

class GenericMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the genericmeasurementcovarianceestimatoroutput operation.
   */
  GenericMeasurementCovarianceEstimatorOutput();
  /**
   * @brief Performs the ~genericmeasurementcovarianceestimatoroutput operation.
   */
  virtual ~GenericMeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the genericmeasurementcovarianceestimatoroutput operation.
   *
   * @param generic_estimator_in The generic estimator.
   */
  GenericMeasurementCovarianceEstimatorOutput(GenericMeasurementCovarianceEstimator* generic_estimator_in);
  
  /**
   * @brief Performs the genericmeasurementcovarianceestimatoroutput operation.
   *
   * @param another The GenericMeasurementCovarianceEstimator instance to copy from.
   */
  GenericMeasurementCovarianceEstimatorOutput(const GenericMeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for GenericMeasurementCovarianceEstimator.
   *
   * @param another The GenericMeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const GenericMeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this GenericMeasurementCovarianceEstimator object.
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
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,const arma::mat &inv_sigma_precomp,
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
  /*
   double evaluate_likelihood(const Parameters &conditioned_on_parameters);
   */
  double subsample_evaluate_likelihood();
  /*
   double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);
   */
  
  MeasurementCovarianceEstimator* get_estimator();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  /** @brief The generic estimator. */
  GenericMeasurementCovarianceEstimator* generic_estimator;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);
  
  //Parameters simulated_measurement;
  /** @brief The measurement state. */
  arma::colvec measurement_state;
  /** @brief The random shift. */
  arma::colvec random_shift;
  
  /**
   * @brief Copies the state of another GenericMeasurementCovarianceEstimator into this object.
   *
   * @param another The GenericMeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const GenericMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
