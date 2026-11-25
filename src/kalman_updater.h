#ifndef KALMANUPDATER_H
#define KALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{

class KalmanFilterOutput;
class KalmanFilter;

class KalmanUpdater
{
  
public:
  
  KalmanUpdater();
  KalmanUpdater(const std::string &state_variable_in,
                const std::string &measurement_variable_in);
  virtual ~KalmanUpdater();
  
  KalmanUpdater(const KalmanUpdater &another);
  
  void operator=(const KalmanUpdater &another);
  virtual KalmanUpdater* duplicate() const=0;
  
  virtual void update(KalmanFilterOutput* current_state,
                      const arma::colvec &current_measurement)=0;
  
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;

  // (C276) the steady-state option of the Kalman filter (see
  // KalmanPredictor): freeze_gain keeps the gain and the factor of the
  // predicted measurement's covariance for the covariance entering this
  // update; update_mean then updates the mean alone with them
  virtual bool supports_steady_state() const { return false; }
  virtual void track_matrix_changes() {}
  virtual bool matrices_changed() const { return true; }
  virtual void freeze_gain(const arma::mat &entering_covariance) { Rcpp::stop("KalmanUpdater::freeze_gain - not available."); }
  virtual void update_mean(KalmanFilterOutput* current_state,
                           const arma::colvec &current_measurement) { Rcpp::stop("KalmanUpdater::update_mean - not available."); }
  
  std::string get_state_variable() const;
  std::string get_measurement_variable() const;
  
protected:
  
  std::string state_variable;
  std::string measurement_variable;
  
  friend KalmanFilter;
  bool set_using_parameters;
  Parameters conditioned_on_parameters;
  
  void make_copy(const KalmanUpdater &another);
  
};
}

#endif
