#ifndef EXACTKALMANPREDICTOR_H
#define EXACTKALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_predictor.h"
#include "ilike_header.h"

namespace ilike
{

class ExactKalmanPredictor : public KalmanPredictor
{
  
public:
  
  ExactKalmanPredictor();
  
  ExactKalmanPredictor(const arma::mat &transition_matrix,
                       const arma::mat &process_noise);
  
  ExactKalmanPredictor(GetMatrixPtr transition_matrix_function,
                       GetMatrixPtr process_noise_function);

  virtual ~ExactKalmanPredictor();
  
  ExactKalmanPredictor(const ExactKalmanPredictor &another);
  
  void operator=(const ExactKalmanPredictor &another);
  KalmanPredictor* duplicate() const;
  
  void predict(KalmanFilterOutput* current_state);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);

  bool supports_steady_state() const { return true; }
  void track_matrix_changes() { this->tracks_matrix_changes = true; }
  bool matrices_changed() const { return this->changed_at_last_set; }
  void predict_mean(KalmanFilterOutput* current_state);
  
protected:
  
  void make_copy(const ExactKalmanPredictor &another);
  
  // Three options:
  // provide matrices
  // provide a function that takes a time step and gives a matrix
  // provide a function that takes parameters and time step and gives a matrix
  
  arma::mat transition_matrix;
  arma::mat process_noise;

  // (C274) workspace, kept between steps so that its memory is reused
  arma::colvec new_mean;
  arma::mat new_covariance;

  // (C276) whether set_parameters compares the new matrices with the last
  // ones, and whether they differed (constant matrices never change)
  bool tracks_matrix_changes = false;
  bool changed_at_last_set = false;
  
  GetMatrixPtr transition_matrix_function;
  GetMatrixPtr process_noise_function;

};
}

#endif
