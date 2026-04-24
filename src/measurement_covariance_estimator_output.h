#ifndef MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file measurement_covariance_estimator_output.h
   * @brief Defines the MeasurementCovarianceEstimator class.
   *
   * Estimates the measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class MeasurementCovarianceEstimator
   * @brief The measurement covariance estimator class.
   */


class MeasurementCovarianceEstimator;
class Index;

class MeasurementCovarianceEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the measurementcovarianceestimatoroutput operation.
   */
  MeasurementCovarianceEstimatorOutput();
  /**
   * @brief Performs the ~measurementcovarianceestimatoroutput operation.
   */
  virtual ~MeasurementCovarianceEstimatorOutput();
  
  /**
   * @brief Performs the measurementcovarianceestimatoroutput operation.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  MeasurementCovarianceEstimatorOutput(const MeasurementCovarianceEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for MeasurementCovarianceEstimator.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const MeasurementCovarianceEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this MeasurementCovarianceEstimator object.
   *
   * @return The result.
   */
  virtual MeasurementCovarianceEstimatorOutput* duplicate() const=0;
  
  /**
   * @brief Simulates the required variables.
   *
   * @param parameters The parameters.
   */
  void simulate(const Parameters &parameters);
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param parameters The parameters.
   */
  void subsample_simulate(const Parameters &parameters);
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param parameters The parameters.
   */
  virtual void specific_simulate(const Parameters &parameters)=0;
  /**
   * @brief Performs the subsample specific simulate operation.
   *
   * @param parameters The parameters.
   */
  virtual void subsample_specific_simulate(const Parameters &parameters)=0;
  /**
   * @brief Returns the measurement state for covariance.
   *
   * @return The result.
   */
  virtual arma::rowvec get_measurement_state_for_covariance() const=0;
  /**
   * @brief Returns the shift.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  virtual arma::colvec get_shift(double inverse_incremental_temperature) const=0;
  /**
   * @brief Returns the deterministic shift.
   *
   * @return The result.
   */
  virtual arma::colvec get_deterministic_shift() const=0;
  
  /**
   * @brief Returns the measurement.
   *
   * @return The result.
   */
  arma::colvec* get_measurement();
  
  virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                    const arma::mat &inv_sigma_precomp,
                                                    double log_det_precomp) const=0;
  /*
   virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                              const arma::mat &inv_sigma_precomp,
                                                              double log_det_precomp) const=0;
  /*
   virtual double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual double evaluate_likelihood()=0;
  //virtual double evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  /**
   * @brief Performs the subsample evaluate likelihood operation.
   *
   * @return The result.
   */
  virtual double subsample_evaluate_likelihood()=0;
  //virtual double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Returns the estimator.
   *
   * @return The result.
   */
  virtual MeasurementCovarianceEstimator* get_estimator()=0;
  
  /**
   * @brief Closes any open file streams.
   */
  virtual void close_ofstreams()=0;
  
  /**
   * @brief Writes results to a file.
   *
   * @param directory_name The directory name.
   */
  void write(const std::string &directory_name);
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index="")=0;
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
protected:
  
  /** @brief The write to file flag. */
  bool write_to_file_flag;
  
  // Stored in ModelAndAlgorithm or in main.
  //MeasurementCovarianceEstimator* estimator;
  
  /**
   * @brief Copies the state of another MeasurementCovarianceEstimator into this object.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const MeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
