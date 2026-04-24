#ifndef DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file direct_linear_gaussian_measurement_covariance_estimator.h
   * @brief Defines the DirectGaussianMeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a direct gaussian measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class DirectGaussianMeasurementCovarianceEstimatorOutput
   * @brief The direct gaussian measurement covariance estimator output class.
   */


//class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectLinearGaussianMeasurementCovarianceEstimator : public DirectGaussianMeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Performs the directlineargaussianmeasurementcovarianceestimator operation.
   */
  DirectLinearGaussianMeasurementCovarianceEstimator();
  
  DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                     size_t* seed_in,
                                                     Data* data_in,
                                                     std::shared_ptr<Transform> transform_in,
                                                     std::shared_ptr<Transform> summary_statistics_in,
                                                     const arma::mat &measurement_matrix_in,
                                                     const arma::mat &measurement_covariance_in,
                                                     const std::string &measurement_variable_in,
                                                     const std::vector<std::string> &state_variables_in);
  
  DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                     size_t* seed_in,
                                                     Data* data_in,
                                                     std::shared_ptr<Transform> transform_in,
                                                     std::shared_ptr<Transform> summary_statistics_in,
                                                     GetMatrixPtr measurement_matrix_in,
                                                     GetMatrixPtr measurement_covariance_in,
                                                     const std::string &measurement_variable_in,
                                                     const std::vector<std::string> &state_variables_in);
  
  /**
   * @brief Performs the ~directlineargaussianmeasurementcovarianceestimator operation.
   */
  virtual ~DirectLinearGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Performs the directlineargaussianmeasurementcovarianceestimator operation.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  DirectLinearGaussianMeasurementCovarianceEstimator(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for DirectGaussianMeasurementCovarianceEstimatorOutput.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  /**
   * @brief Creates a deep copy of this DirectGaussianMeasurementCovarianceEstimatorOutput object.
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
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  /**
   * @brief Returns the measurement state parameters.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  Parameters get_measurement_state_parameters(const Parameters &parameters) const;
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Dhathalf,
                           const arma::mat &P,
                           const arma::mat &Vtranspose,
                           const arma::mat &Yhat,
                           double inverse_incremental_temperature) const;
  
  arma::mat get_sqrt_adjustment(const arma::mat &Sigma,
                                const arma::mat &HSigmaHt,
                                double inverse_incremental_temperature) const;
  
  //GaussianIndependentProposalKernel get_kernel() const;
  
  //arma::mat get_Cygivenx() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  //void setup();
  //void setup(const Parameters &parameters);
  
protected:
  
  //friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
  //MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  //MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
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
   * @brief Copies the state of another DirectGaussianMeasurementCovarianceEstimatorOutput into this object.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void make_copy(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  /** @brief The as. */
  std::vector<arma::mat> As;
  /** @brief The a functions. */
  std::vector<GetMatrixPtr> A_functions;
  /** @brief The measurement noise functions. */
  std::vector<GetMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;
  
  /** @brief The state variables. */
  std::vector<std::string> state_variables;
  
};
}

#endif
