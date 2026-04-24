#ifndef EXACTKALMANPREDICTOR_H
#define EXACTKALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_predictor.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file exact_kalman_predictor.h
   * @brief Defines the ExactKalmanPredictor class.
   *
   * Implements exact kalman predictor. Performs Gaussian state-space inference using Kalman recursions.
   *
   * @namespace ilike
   * @class ExactKalmanPredictor
   * @brief An exact kalman predictor derived from KalmanPredictor.
   */


class ExactKalmanPredictor : public KalmanPredictor
{
  
public:
  
  /**
   * @brief Default constructor for ExactKalmanPredictor.
   */
  ExactKalmanPredictor();
  
  ExactKalmanPredictor(const arma::mat &transition_matrix,
                       const arma::mat &process_noise);
  
  ExactKalmanPredictor(GetMatrixPtr transition_matrix_function,
                       GetMatrixPtr process_noise_function);
  
  //ExactKalmanPredictor(GetProcessMatrixFromParametersPtr transition_matrix_parameters_function,
  //                     GetProcessMatrixFromParametersPtr process_noise_parameters_function);
  
  //ExactKalmanPredictor(GetProcessMatrixFromTimeParametersPtr transition_matrix_time_parameters_function,
  //                     GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function);
  
  /**
   * @brief Destructor for ExactKalmanPredictor.
   */
  virtual ~ExactKalmanPredictor();
  
  /**
   * @brief Copy constructor for ExactKalmanPredictor.
   *
   * @param another The ExactKalmanPredictor instance to copy from.
   */
  ExactKalmanPredictor(const ExactKalmanPredictor &another);
  
  /**
   * @brief Assignment operator for ExactKalmanPredictor.
   *
   * @param another The ExactKalmanPredictor instance to copy from.
   */
  void operator=(const ExactKalmanPredictor &another);
  /**
   * @brief Creates a deep copy of this ExactKalmanPredictor object.
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
   * @brief Copies the state of another ExactKalmanPredictor into this object.
   *
   * @param another The ExactKalmanPredictor instance to copy from.
   */
  void make_copy(const ExactKalmanPredictor &another);
  
  // Three options:
  // provide matrices
  // provide a function that takes a time step and gives a matrix
  // provide a function that takes parameters and time step and gives a matrix
  
  /** @brief The transition matrix. */
  arma::mat transition_matrix;
  /** @brief The process noise. */
  arma::mat process_noise;
  
  /** @brief The transition matrix function. */
  GetMatrixPtr transition_matrix_function;
  /** @brief The process noise function. */
  GetMatrixPtr process_noise_function;
  
  //GetProcessMatrixFromTimePtr transition_matrix_time_function;
  //GetProcessMatrixFromTimePtr process_noise_time_function;
  
  //GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
  //GetProcessMatrixFromParametersPtr process_noise_parameters_function;
  
  //GetProcessMatrixFromTimeParametersPtr transition_matrix_time_parameters_function;
  //GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function;
  
};
}

#endif
