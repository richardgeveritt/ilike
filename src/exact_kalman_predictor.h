#ifndef EXACTKALMANPREDICTOR_H
#define EXACTKALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_predictor.h"
#include "ilike_header.h"

class ExactKalmanPredictor : public KalmanPredictor
{

public:

  ExactKalmanPredictor();
  
  ExactKalmanPredictor(const arma::mat &transition_matrix,
                       const arma::mat &process_noise);
  
  ExactKalmanPredictor(GetProcessMatrixFromTimePtr transition_matrix_time_function,
                       GetProcessMatrixFromTimePtr process_noise_time_function);
  
  //ExactKalmanPredictor(GetProcessMatrixFromParametersPtr transition_matrix_parameters_function,
  //                     GetProcessMatrixFromParametersPtr process_noise_parameters_function);
  
  ExactKalmanPredictor(GetProcessMatrixFromTimeParametersPtr transition_matrix_time_parameters_function,
                       GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function);

  virtual ~ExactKalmanPredictor();

  ExactKalmanPredictor(const ExactKalmanPredictor &another);

  void operator=(const ExactKalmanPredictor &another);
  KalmanPredictor* duplicate() const;

  void predict(KalmanFilterOutput* current_state,
               double current_time,
               double next_time);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
protected:
  
  void make_copy(const ExactKalmanPredictor &another);
  
  // Three options:
  // provide matrices
  // provide a function that takes a time step and gives a matrix
  // provide a function that takes parameters and time step and gives a matrix
  
  arma::mat transition_matrix;
  arma::mat process_noise;
  
  GetProcessMatrixFromTimePtr transition_matrix_time_function;
  GetProcessMatrixFromTimePtr process_noise_time_function;
  
  //GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
  //GetProcessMatrixFromParametersPtr process_noise_parameters_function;
  
  GetProcessMatrixFromTimeParametersPtr transition_matrix_time_parameters_function;
  GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function;

};

#endif
