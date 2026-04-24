#ifndef UNSCENTEDKALMANUPDATER_H
#define UNSCENTEDKALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_updater.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file unscented_kalman_updater.h
   * @brief Defines the KalmanFilterOutput class.
   *
   * Stores and manages the output produced by KalmanFilter. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class KalmanFilterOutput
   * @brief The kalman filter output class.
   */


class KalmanFilterOutput;

class UnscentedKalmanUpdater : public KalmanUpdater
{
  
public:
  
  /**
   * @brief Performs the unscentedkalmanupdater operation.
   */
  UnscentedKalmanUpdater();
  
  UnscentedKalmanUpdater(const std::string &state_variable_in,
                         const std::string &measurement_variable_in,
                         GetMatrixSimulateMeasurementKernelPtr measurement_kernel_function_in,
                         GetMatrixPtr measurement_noise_function_in,
                         double w0_in = 1.0/3.0);
  
  UnscentedKalmanUpdater(const std::string &state_variable_in,
                         const std::string &measurement_variable_in,
                         MatrixSimulateMeasurementKernelPtr measurement_kernel_in,
                         const arma::mat &measurement_noise_in,
                         double w0_in = 1.0/3.0);
  
  /**
   * @brief Performs the ~unscentedkalmanupdater operation.
   */
  virtual ~UnscentedKalmanUpdater();
  
  /**
   * @brief Performs the unscentedkalmanupdater operation.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  UnscentedKalmanUpdater(const UnscentedKalmanUpdater &another);
  
  /**
   * @brief Assignment operator for KalmanFilterOutput.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void operator=(const UnscentedKalmanUpdater &another);
  /**
   * @brief Creates a deep copy of this KalmanFilterOutput object.
   *
   * @return The result.
   */
  KalmanUpdater* duplicate() const;
  
  void update(KalmanFilterOutput* current_state,
              const arma::colvec &current_measurement);
  
  /**
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
protected:
  
  /**
   * @brief Copies the state of another KalmanFilterOutput into this object.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void make_copy(const UnscentedKalmanUpdater &another);
  
  /** @brief The w0. */
  double w0;
  
  /** @brief The conditioned on parameters. */
  Parameters conditioned_on_parameters;
  
  // use either these
  /** @brief The measurement kernel function. */
  GetMatrixSimulateMeasurementKernelPtr measurement_kernel_function;
  /** @brief The measurement noise function. */
  GetMatrixPtr measurement_noise_function;
  
  // or these
  /** @brief The measurement kernel. */
  MatrixSimulateMeasurementKernelPtr measurement_kernel;
  /** @brief The measurement noise. */
  arma::mat measurement_noise;
  
};
}

#endif
