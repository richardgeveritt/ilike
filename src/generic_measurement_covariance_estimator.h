#ifndef GENERICMEASUREMENTCOVARIANCEESTIMATOR_H
#define GENERICMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file generic_measurement_covariance_estimator.h
   * @brief Defines the MeasurementCovarianceEstimatorOutput class.
   *
   * Stores the output of a measurement covariance estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class MeasurementCovarianceEstimatorOutput
   * @brief The measurement covariance estimator output class.
   */


class MeasurementCovarianceEstimatorOutput;
class GenericMeasurementCovarianceEstimatorOutput;

class GenericMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{
  
public:
  
  /**
   * @brief Performs the genericmeasurementcovarianceestimator operation.
   */
  GenericMeasurementCovarianceEstimator();
  
  GenericMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                        size_t* seed_in,
                                        Data* data_in,
                                        std::shared_ptr<Transform> transform_in,
                                        std::shared_ptr<Transform> summary_statistics_in,
                                        SimulateModelPtr simulator_in,
                                        const std::vector<std::string> &measurement_variables_in);
  
  /**
   * @brief Performs the ~genericmeasurementcovarianceestimator operation.
   */
  virtual ~GenericMeasurementCovarianceEstimator();
  
  /**
   * @brief Performs the genericmeasurementcovarianceestimator operation.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  GenericMeasurementCovarianceEstimator(const GenericMeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for MeasurementCovarianceEstimatorOutput.
   *
   * @param another The MeasurementCovarianceEstimatorOutput instance to copy from.
   */
  void operator=(const GenericMeasurementCovarianceEstimator &another);
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
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  /**
   * @brief Performs the need cxx operation.
   *
   * @return The result.
   */
  bool need_Cxx() const;
  
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
  
  friend GenericMeasurementCovarianceEstimatorOutput;
  
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
  void make_copy(const GenericMeasurementCovarianceEstimator &another);
  
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
  
  /** @brief The dimension. */
  size_t dimension;
  
  /** @brief The cygivenx. */
  arma::mat Cygivenx;
  
};
}

#endif
