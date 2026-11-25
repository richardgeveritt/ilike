#ifndef EXACTKALMANUPDATER_H
#define EXACTKALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_updater.h"
#include "ilike_header.h"

namespace ilike
{

class ExactKalmanUpdater : public KalmanUpdater
{
  
public:
  
  ExactKalmanUpdater();
  
  ExactKalmanUpdater(const std::string &state_variable_in,
                     const std::string &measurement_variable_in,
                     const arma::mat &measurement_matrix_in,
                     const arma::mat &measurement_noise_in);
  
  ExactKalmanUpdater(const std::string &state_variable_in,
                     const std::string &measurement_variable_in,
                     GetMatrixPtr measurement_matrix_function,
                     GetMatrixPtr measurement_noise_function);
  
  virtual ~ExactKalmanUpdater();
  
  ExactKalmanUpdater(const ExactKalmanUpdater &another);
  
  void operator=(const ExactKalmanUpdater &another);
  KalmanUpdater* duplicate() const;
  
  void update(KalmanFilterOutput* current_state,
              const arma::colvec &current_measurement);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);

  bool supports_steady_state() const { return true; }
  void track_matrix_changes() { this->tracks_matrix_changes = true; }
  bool matrices_changed() const { return this->changed_at_last_set; }
  void freeze_gain(const arma::mat &entering_covariance);
  void update_mean(KalmanFilterOutput* current_state,
                   const arma::colvec &current_measurement);
  
protected:
  
  void make_copy(const ExactKalmanUpdater &another);
  
  GetMatrixPtr measurement_matrix_function;
  GetMatrixPtr measurement_noise_function;
  
  arma::mat measurement_matrix;
  arma::mat measurement_noise;

  // (C274) workspace, kept between steps so that its memory is reused (the
  // values are not state: each update overwrites them)
  arma::colvec innovation;
  arma::mat G, S, L, W;
  arma::colvec z;
  arma::colvec new_mean;
  arma::mat new_covariance;
  arma::rowvec g;

  // (C276) as in ExactKalmanPredictor; and the frozen gain K, the factor L
  // of the frozen predicted measurement's covariance, and the constant of
  // its log-density, -(d/2) log(2 pi) - sum(log diag L)
  bool tracks_matrix_changes = false;
  bool changed_at_last_set = false;
  arma::mat frozen_gain;
  arma::mat frozen_L;
  double frozen_log_constant = 0.0;
  // (C275) the components of a measurement with a diagonal noise covariance
  // taken one at a time (scalar updates; exact)
  double update_one_component_at_a_time(const arma::colvec &predicted_mean,
                                        const arma::mat &predicted_covariance,
                                        const arma::colvec &current_measurement);
};
}

#endif
