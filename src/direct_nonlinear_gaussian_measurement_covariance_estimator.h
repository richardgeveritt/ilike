#ifndef DIRECTNONLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTNONLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file direct_nonlinear_gaussian_measurement_covariance_estimator.h
   * @brief Defines the DirectGaussianMeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a direct gaussian measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class DirectGaussianMeasurementCovarianceEstimatorOutput
   * @brief The direct gaussian measurement covariance estimator output class.
   */


//class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectNonLinearGaussianMeasurementCovarianceEstimator : public DirectGaussianMeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Performs the directnonlineargaussianmeasurementcovarianceestimator operation.
   */
  DirectNonLinearGaussianMeasurementCovarianceEstimator();
  
  
  /*
   DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   std::shared_ptr<Transform> transform_in,
   std::shared_ptr<Transform> summary_statistics_in,
   std::shared_ptr<Transform> transform_function_in,
   const arma::mat &measurement_noise_in);
   //const std::string &measurement_variable_in);
   
   DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   std::shared_ptr<Transform> transform_in,
   std::shared_ptr<Transform> summary_statistics_in,
   std::shared_ptr<Transform> transform_function_in,
   const std::vector<arma::mat> &measurement_noise_in);
   //const std::vector<std::string> &measurement_variable_in);
   
   DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   std::shared_ptr<Transform> transform_in,
   std::shared_ptr<Transform> summary_statistics_in,
   std::shared_ptr<Transform> transform_function_in,
   GetMatrixPtr measurement_noise_function_in);
   //const std::string &measurement_variable_in);
   
   DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   std::shared_ptr<Transform> transform_in,
   std::shared_ptr<Transform> summary_statistics_in,
   std::shared_ptr<Transform> transform_function_in,
   const std::vector<GetMatrixPtr> &measurement_noise_functions_in);
   //const std::vector<std::string> &measurements_variable_in);
   */
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const arma::mat &measurement_noise_in,
                                                        const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<arma::mat> &measurement_noise_in,
                                                        const std::vector<std::string> &measurement_variables_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        GetMatrixPtr measurement_noise_function_in,
                                                        const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<GetMatrixPtr> &measurement_noise_functions_in,
                                                        const std::vector<std::string> &measurements_variables_in);
  
  /**
   * @brief Performs the ~directnonlineargaussianmeasurementcovarianceestimator operation.
   */
  virtual ~DirectNonLinearGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Performs the directnonlineargaussianmeasurementcovarianceestimator operation.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  DirectNonLinearGaussianMeasurementCovarianceEstimator(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for DirectGaussianMeasurementCovarianceEstimatorOutput.
   *
   * @param another The DirectGaussianMeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);
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
  void make_copy(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  /** @brief The transform function. */
  std::shared_ptr<Transform> transform_function;
  /** @brief The measurement noise functions. */
  std::vector<GetMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;
  
};
}

#endif
