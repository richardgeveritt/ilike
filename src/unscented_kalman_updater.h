#ifndef UNSCENTEDKALMANUPDATER_H
#define UNSCENTEDKALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "kalman_updater.h"
#include "ilike_header.h"

namespace ilike
{
class KalmanFilterOutput;

class UnscentedKalmanUpdater : public KalmanUpdater
{
  
public:
  
  UnscentedKalmanUpdater();
  
  UnscentedKalmanUpdater(const std::string &state_variable_in,
                         const std::string &measurement_variable_in,
                         GetMatrixSimulateMeasurementKernelPtr measurement_kernel_function_in,
                         GetMatrixPtr measurement_noise_function_in,
                         double w0_in = 1.0/3.0);
  
  UnscentedKalmanUpdater(const std::string &state_variable_in,
                         const std::string &measurement_variable_in,
                         MatrixSimulateMeasurementKernelPtr measurement_kernel_in,
                         const arma::mat &measurement_noise_in,
                         double w0_in = 1.0/3.0);
  
  virtual ~UnscentedKalmanUpdater();
  
  UnscentedKalmanUpdater(const UnscentedKalmanUpdater &another);
  
  void operator=(const UnscentedKalmanUpdater &another);
  KalmanUpdater* duplicate() const;
  
  void update(KalmanFilterOutput* current_state,
              const arma::colvec &current_measurement);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
protected:
  
  void make_copy(const UnscentedKalmanUpdater &another);
  
  double w0;
  
  Parameters conditioned_on_parameters;
  
  // use either these
  GetMatrixSimulateMeasurementKernelPtr measurement_kernel_function;
  GetMatrixPtr measurement_noise_function;
  
  // or these
  MatrixSimulateMeasurementKernelPtr measurement_kernel;
  arma::mat measurement_noise;
  
};
}

#endif
