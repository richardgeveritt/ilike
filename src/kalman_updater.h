#ifndef KALMANUPDATER_H
#define KALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{
  /**
   * @file kalman_updater.h
   * @brief Defines the KalmanFilterOutput class.
   *
   * Stores and manages the output produced by KalmanFilter. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class KalmanFilterOutput
   * @brief The kalman filter output class.
   */


class KalmanFilterOutput;
class KalmanFilter;

class KalmanUpdater
{
  
public:
  
  /**
   * @brief Performs the kalmanupdater operation.
   */
  KalmanUpdater();
  KalmanUpdater(const std::string &state_variable_in,
                const std::string &measurement_variable_in);
  /**
   * @brief Performs the ~kalmanupdater operation.
   */
  virtual ~KalmanUpdater();
  
  /**
   * @brief Performs the kalmanupdater operation.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  KalmanUpdater(const KalmanUpdater &another);
  
  /**
   * @brief Assignment operator for KalmanFilterOutput.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void operator=(const KalmanUpdater &another);
  /**
   * @brief Creates a deep copy of this KalmanFilterOutput object.
   *
   * @return The result.
   */
  virtual KalmanUpdater* duplicate() const=0;
  
  virtual void update(KalmanFilterOutput* current_state,
                      const arma::colvec &current_measurement)=0;
  
  /**
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  /**
   * @brief Returns the state variable.
   *
   * @return The result.
   */
  std::string get_state_variable() const;
  /**
   * @brief Returns the measurement variable.
   *
   * @return The result.
   */
  std::string get_measurement_variable() const;
  
protected:
  
  /** @brief The state variable. */
  std::string state_variable;
  /** @brief The measurement variable. */
  std::string measurement_variable;
  
  friend KalmanFilter;
  /** @brief The set using parameters. */
  bool set_using_parameters;
  /** @brief The conditioned on parameters. */
  Parameters conditioned_on_parameters;
  
  /**
   * @brief Copies the state of another KalmanFilterOutput into this object.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void make_copy(const KalmanUpdater &another);
  
};
}

#endif
