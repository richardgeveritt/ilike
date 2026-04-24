#ifndef MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file mixed_generic_direct_gaussian_measurement_covariance_estimator.h
   * @brief Defines the MeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class MeasurementCovarianceEstimatorOutput
   * @brief The measurement covariance estimator output class.
   */


class MeasurementCovarianceEstimatorOutput;
class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;

class MixedGenericDirectGaussianMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Performs the mixedgenericdirectgaussianmeasurementcovarianceestimator operation.
   */
  MixedGenericDirectGaussianMeasurementCovarianceEstimator();
  
  MixedGenericDirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                           size_t* seed_in,
                                                           Data* data_in,
                                                           std::shared_ptr<Transform> transform_in,
                                                           std::shared_ptr<Transform> summary_statistics_in,
                                                           Data* prior_data_in,
                                                           SimulateModelPtr simulator_in,
                                                           const std::vector<std::string> &prior_measurement_variables_in,
                                                           const std::vector<arma::mat> &prior_measurement_noises_in,
                                                           const std::vector<std::string> &data_measurement_variables_in);
  
  /**
   * @brief Performs the ~mixedgenericdirectgaussianmeasurementcovarianceestimator operation.
   */
  virtual ~MixedGenericDirectGaussianMeasurementCovarianceEstimator();
  
  /**
   * @brief Performs the mixedgenericdirectgaussianmeasurementcovarianceestimator operation.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  MixedGenericDirectGaussianMeasurementCovarianceEstimator(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for MeasurementCovarianceEstimatorOutput.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);
  /**
   * @brief Creates a deep copy of this MeasurementCovarianceEstimatorOutput object.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimator* duplicate() const;
  
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
  
  /**
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  /**
   * @brief Performs the need cxx operation.
   *
   * @return The result.
   */
  bool need_Cxx() const;
  
  void find_partial_Cygivenx(const arma::mat &Cxy,
                             const arma::mat &Cyy,
                             const arma::mat &packed_members);
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
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
   * @brief Returns the cygivenx.
   *
   * @return The result.
   */
  arma::mat get_Cygivenx() const;
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature) const;
  
  /**
   * @brief Simulates the required variables.
   *
   * @param current_state The current state.
   *
   * @return The result.
   */
  Parameters simulate(const Parameters &current_state);
  
  /**
   * @brief Performs the gaussian simulate operation.
   *
   * @return The result.
   */
  arma::colvec gaussian_simulate();
  
  /**
   * @brief Returns the prior measurement covariance.
   *
   * @return The result.
   */
  arma::mat get_prior_measurement_covariance() const;
  /**
   * @brief Returns the prior measurement covariance embedded in full space.
   *
   * @return The result.
   */
  arma::mat get_prior_measurement_covariance_embedded_in_full_space() const;
  /**
   * @brief Returns the measurement covariance for likelihood ratio.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  arma::mat get_measurement_covariance_for_likelihood_ratio(double inverse_incremental_temperature) const;
  
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
  
  friend MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
  
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
   * @brief Copies the state of another MeasurementCovarianceEstimatorOutput into this object.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  //GetMeasurementMatrixPtr measurement_noise_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //arma::mat measurement_noise;
  
  /** @brief The simulator. */
  SimulateModelPtr simulator;
  
  /** @brief The gaussian simulator. */
  GaussianIndependentProposalKernel gaussian_simulator;
  
  /** @brief The measurement dimension. */
  size_t measurement_dimension; // dimension of prior
  
  /** @brief The cygivenx. */
  arma::mat Cygivenx;
  /** @brief The likelihood cygivenx. */
  arma::mat likelihood_Cygivenx;
  
  // mean of zero
  /** @brief The kernel. */
  GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  /** @brief The transform function. */
  TransformToDataPtr transform_function;
  /** @brief The prior measurement noise functions. */
  std::vector<GetMatrixPtr> prior_measurement_noise_functions;
  
  // in this class, this->measurement_variables (in base) always refers to the "true" measurement variables, for the likelihood term - not those for the prior, which we store in this derived class (same goes for the prior "Data")
  /** @brief The prior measurement variables. */
  std::vector<std::string> prior_measurement_variables;
  
  // not stored here - points to prior state mean in base
  /** @brief The prior data. */
  Data* prior_data;
  
};
}

#endif
