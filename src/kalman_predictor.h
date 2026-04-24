#ifndef KALMANPREDICTOR_H
#define KALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{
  /**
   * @file kalman_predictor.h
   * @brief Defines the KalmanFilterOutput class.
   *
   * Stores and manages the output produced by KalmanFilter. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class KalmanFilterOutput
   * @brief The kalman filter output class.
   */


class KalmanFilterOutput;

class KalmanPredictor
{
  
public:
  
  /**
   * @brief Performs the kalmanpredictor operation.
   */
  KalmanPredictor();
  /**
   * @brief Performs the ~kalmanpredictor operation.
   */
  virtual ~KalmanPredictor();
  
  /**
   * @brief Performs the kalmanpredictor operation.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  KalmanPredictor(const KalmanPredictor &another);
  
  /**
   * @brief Assignment operator for KalmanFilterOutput.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void operator=(const KalmanPredictor &another);
  /**
   * @brief Creates a deep copy of this KalmanFilterOutput object.
   *
   * @return The result.
   */
  virtual KalmanPredictor* duplicate() const=0;
  
  /**
   * @brief Performs the predict operation.
   *
   * @param current_state The current state.
   */
  virtual void predict(KalmanFilterOutput* current_state)=0;
  
  //virtual void predict(KalmanFilterOutput* current_state,
  //                     double current_time,
  //                     double next_time)=0;
  
  /**
   * @brief Sets the parameters.
   *
   * @param conditioned_on_parameters_in The conditioned on parameters.
   */
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
protected:
  
  /** @brief The set using parameters. */
  bool set_using_parameters;
  //bool set_using_time;
  /** @brief The conditioned on parameters. */
  Parameters conditioned_on_parameters;
  
  /**
   * @brief Copies the state of another KalmanFilterOutput into this object.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void make_copy(const KalmanPredictor &another);
  
};
}

#endif
