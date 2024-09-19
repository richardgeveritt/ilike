#ifndef UNSCENTEDKALMANPREDICTOR_H
#define UNSCENTEDKALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_predictor.h"
#include "ilike_header.h"

namespace ilike
{
class KalmanFilterOutput;

class UnscentedKalmanPredictor : public KalmanPredictor
{
  
public:
  
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
  
  virtual ~UnscentedKalmanPredictor();
  
  UnscentedKalmanPredictor(const UnscentedKalmanPredictor &another);
  
  void operator=(const UnscentedKalmanPredictor &another);
  KalmanPredictor* duplicate() const;
  
  void predict(KalmanFilterOutput* current_state);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
protected:
  
  void make_copy(const UnscentedKalmanPredictor &another);
  
  double w0;
  
  // Three options:
  // provide kernel/matrix
  // provide a function that takes a time step and gives a kernel/matrix
  // provide a function that takes parameters and time step and gives a kernel/matrix
  
  SimulateTransitionKernelPtr transition_kernel;
  arma::mat process_noise;
  
  SimulateTransitionKernelFromParametersPtr transition_kernel_function;
  GetMatrixPtr process_noise_function;
  
  //GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
  //GetProcessMatrixFromParametersPtr process_noise_parameters_function;
  
  //SimulateTransitionKernelFromTimeParametersPtr transition_kernel_time_parameters_function;
  //GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function;
  
};
}

#endif
