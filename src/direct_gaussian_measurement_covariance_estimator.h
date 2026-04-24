#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file direct_gaussian_measurement_covariance_estimator.h
   * @brief Defines the DirectGaussianMeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a direct gaussian measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class DirectGaussianMeasurementCovarianceEstimatorOutput
   * @brief The direct gaussian measurement covariance estimator output class.
   */


class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectGaussianMeasurementCovarianceEstimator : public GaussianMeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Performs the directgaussianmeasurementcovarianceestimator operation.
   */
  DirectGaussianMeasurementCovarianceEstimator();
  
  /*
   DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   std::shared_ptr<Transform> transform_in,
   std::shared_ptr<Transform> summary_statistics_in,
   const std::string &measurement_variable_in);
   */
  
  DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               std::shared_ptr<Transform> transform_in,
                                               std::shared_ptr<Transform> summary_statistics_in,
                                               const std::vector<std::string> &measurement_variables_in);
  
  /**
   * @brief Performs the ~directgaussianmeasurementcovarianceestimator operation.
   */
  virtual ~DirectGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Performs the directgaussianmeasurementcovarianceestimator operation.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for DirectGaussianMeasurementCovarianceEstimatorOutput.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  /**
   * @brief Returns the measurement state.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  arma::colvec get_measurement_state(const Parameters &parameters) const;
  /**
   * @brief Returns the measurement state parameters.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  virtual Parameters get_measurement_state_parameters(const Parameters &parameters) const=0;
  
  //virtual GaussianIndependentProposalKernel get_kernel() const=0;
  
  /**
   * @brief Returns the cygivenx.
   *
   * @return The result.
   */
  arma::mat get_Cygivenx() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  /**
   * @brief Returns the measurement covariance.
   *
   * @return The result.
   */
  arma::mat get_measurement_covariance() const;
  
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
  
  // mean of zero
  /** @brief The kernel. */
  GaussianIndependentProposalKernel kernel;
  
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
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
  
  //void setup_measurement_variables();
  //void setup_measurement_variables(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Copies the state of another DirectGaussianMeasurementCovarianceEstimatorOutput into this object.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void make_copy(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //std::shared_ptr<Transform> transform_function;
  //std::vector<GetMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;
  
};
}

#endif
