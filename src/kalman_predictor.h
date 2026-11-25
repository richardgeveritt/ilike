#ifndef KALMANPREDICTOR_H
#define KALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{

class KalmanFilterOutput;

class KalmanPredictor
{
  
public:
  
  KalmanPredictor();
  virtual ~KalmanPredictor();
  
  KalmanPredictor(const KalmanPredictor &another);
  
  void operator=(const KalmanPredictor &another);
  virtual KalmanPredictor* duplicate() const=0;
  
  virtual void predict(KalmanFilterOutput* current_state)=0;

  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;

  // (C276) the steady-state option of the Kalman filter: whether this
  // predictor takes part (the default: no, and the filter never freezes);
  // whether its matrices changed at the last set_parameters (tracked once
  // track_matrix_changes() has been called); and the prediction of the mean
  // alone (the covariance is frozen)
  virtual bool supports_steady_state() const { return false; }
  virtual void track_matrix_changes() {}
  virtual bool matrices_changed() const { return true; }
  virtual void predict_mean(KalmanFilterOutput* current_state) { Rcpp::stop("KalmanPredictor::predict_mean - not available."); }

protected:
  
  bool set_using_parameters;
  Parameters conditioned_on_parameters;
  
  void make_copy(const KalmanPredictor &another);
  
};
}

#endif
