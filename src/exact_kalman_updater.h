#ifndef EXACTKALMANUPDATER_H
#define EXACTKALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_updater.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file exact_kalman_updater.h
   * @brief Defines the ExactKalmanUpdater class.
   *
   * Implements exact kalman updater. Performs Gaussian state-space inference using Kalman recursions.
   *
   * @namespace ilike
   * @class ExactKalmanUpdater
   * @brief An exact kalman updater derived from KalmanUpdater.
   */


class ExactKalmanUpdater : public KalmanUpdater
{
  
public:
  
  /**
   * @brief Default constructor for ExactKalmanUpdater.
   */
  ExactKalmanUpdater();
  
  ExactKalmanUpdater(const std::string &state_variable_in,
                     const std::string &measurement_variable_in,
                     const arma::mat &measurement_matrix_in,
                     const arma::mat &measurement_noise_in);
  
  ExactKalmanUpdater(const std::string &state_variable_in,
                     const std::string &measurement_variable_in,
                     GetMatrixPtr measurement_matrix_function,
                     GetMatrixPtr measurement_noise_function);
  
  /**
   * @brief Destructor for ExactKalmanUpdater.
   */
  virtual ~ExactKalmanUpdater();
  
  /**
   * @brief Copy constructor for ExactKalmanUpdater.
   *
   * @param another The ExactKalmanUpdater instance to copy from.
   */
  ExactKalmanUpdater(const ExactKalmanUpdater &another);
  
  /**
   * @brief Assignment operator for ExactKalmanUpdater.
   *
   * @param another The ExactKalmanUpdater instance to copy from.
   */
  void operator=(const ExactKalmanUpdater &another);
  /**
   * @brief Creates a deep copy of this ExactKalmanUpdater object.
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
   * @brief Copies the state of another ExactKalmanUpdater into this object.
   *
   * @param another The ExactKalmanUpdater instance to copy from.
   */
  void make_copy(const ExactKalmanUpdater &another);
  
  /** @brief The measurement matrix function. */
  GetMatrixPtr measurement_matrix_function;
  /** @brief The measurement noise function. */
  GetMatrixPtr measurement_noise_function;
  
  /** @brief The measurement matrix. */
  arma::mat measurement_matrix;
  /** @brief The measurement noise. */
  arma::mat measurement_noise;
  
};
}

#endif
