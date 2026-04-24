#ifndef UNSCENTEDKALMANPREDICTOR_H
#define UNSCENTEDKALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_predictor.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file unscented_kalman_predictor.h
   * @brief Defines the KalmanFilterOutput class.
   *
   * Stores and manages the output produced by KalmanFilter. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class KalmanFilterOutput
   * @brief The kalman filter output class.
   */


class KalmanFilterOutput;

class UnscentedKalmanPredictor : public KalmanPredictor
{
  
public:
  
  /**
   * @brief Performs the unscentedkalmanpredictor operation.
   */
  UnscentedKalmanPredictor();
  
  UnscentedKalmanPredictor(SimulateTransitionKernelPtr transition_kernel_in,
                           const arma::mat &process_noise_in,
                           double w0_in = 1.0/3.0);
  
  UnscentedKalmanPredictor(SimulateTransitionKernelFromParametersPtr transition_kernel_function_in,
                           GetMatrixPtr process_noise_function_in,
                           double w0_in = 1.0/3.0);
  
  //GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
  //GetProcessMatrixFromParametersPtr process_noise_parameters_function;
  
  //UnscentedKalmanPredictor(SimulateTransitionKernelFromTimeParametersPtr transition_kernel_time_parameters_function_in,
  //                         GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function_in,
  //                         double w0_in = 1.0/3.0);
  
  /**
   * @brief Performs the ~unscentedkalmanpredictor operation.
   */
  virtual ~UnscentedKalmanPredictor();
  
  /**
   * @brief Performs the unscentedkalmanpredictor operation.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  UnscentedKalmanPredictor(const UnscentedKalmanPredictor &another);
  
  /**
   * @brief Assignment operator for KalmanFilterOutput.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void operator=(const UnscentedKalmanPredictor &another);
  /**
   * @brief Creates a deep copy of this KalmanFilterOutput object.
   *
   * @return The result.
   */
  KalmanPredictor* duplicate() const;
  
  /**
   * @brief Performs the predict operation.
   *
   * @param current_state The current state.
   */
  void predict(KalmanFilterOutput* current_state);
  
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
  void make_copy(const UnscentedKalmanPredictor &another);
  
  /** @brief The w0. */
  double w0;
  
  // Three options:
  // provide kernel/matrix
  // provide a function that takes a time step and gives a kernel/matrix
  // provide a function that takes parameters and time step and gives a kernel/matrix
  
  /** @brief The transition kernel. */
  SimulateTransitionKernelPtr transition_kernel;
  /** @brief The process noise. */
  arma::mat process_noise;
  
  /** @brief The transition kernel function. */
  SimulateTransitionKernelFromParametersPtr transition_kernel_function;
  /** @brief The process noise function. */
  GetMatrixPtr process_noise_function;
  
  //GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
  //GetProcessMatrixFromParametersPtr process_noise_parameters_function;
  
  //SimulateTransitionKernelFromTimeParametersPtr transition_kernel_time_parameters_function;
  //GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function;
  
};
}

#endif
