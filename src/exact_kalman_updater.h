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
  
protected:
  
  void make_copy(const ExactKalmanUpdater &another);
  
  GetMatrixPtr measurement_matrix_function;
  GetMatrixPtr measurement_noise_function;
  
  arma::mat measurement_matrix;
  arma::mat measurement_noise;
  
};
}

#endif
