#include "exact_kalman_updater.h"
#include "utils.h"
#include "kalman_filter_output.h"

ExactKalmanUpdater::ExactKalmanUpdater()
  :KalmanUpdater()
{
}

ExactKalmanUpdater::ExactKalmanUpdater(const arma::mat &measurement_matrix_in,
                                       const arma::mat &measurement_noise_in)
{
  this->measurement_matrix = measurement_matrix_in;
  this->measurement_noise = measurement_noise_in;
  this->set_using_parameters = false;
}

ExactKalmanUpdater::ExactKalmanUpdater(GetMeasurementMatrixPtr measurement_matrix_function_in,
                                       GetMeasurementMatrixPtr measurement_noise_function_in)
{
  this->measurement_matrix_function = measurement_matrix_function_in;
  this->measurement_noise_function = measurement_noise_function_in;
  this->set_using_parameters = true;
}

ExactKalmanUpdater::~ExactKalmanUpdater()
{
  
}

ExactKalmanUpdater::ExactKalmanUpdater(const ExactKalmanUpdater &another)
  :KalmanUpdater(another)
{
  this->make_copy(another);
}

void ExactKalmanUpdater::operator=(const ExactKalmanUpdater &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  KalmanUpdater::operator=(another);
  this->make_copy(another);
}

KalmanUpdater* ExactKalmanUpdater::duplicate() const
{
  return( new ExactKalmanUpdater(*this));
}

void ExactKalmanUpdater::make_copy(const ExactKalmanUpdater &another)
{
  this->measurement_matrix = another.measurement_matrix;
  this->measurement_noise = another.measurement_noise;
  this->measurement_matrix_function = another.measurement_matrix_function;
  this->measurement_noise_function = another.measurement_noise_function;
}

void ExactKalmanUpdater::update(KalmanFilterOutput* current_state,
                                const arma::colvec &current_measurement)
{
  arma::colvec predicted_mean = current_state->predicted_mean_back();
  arma::mat predicted_covariance = current_state->predicted_covariance_back();
  arma::colvec predicted_measurement = this->measurement_matrix*predicted_mean;
  arma::colvec innovation = current_measurement - predicted_measurement;
  arma::mat innovation_covariance = this->measurement_matrix*predicted_covariance*this->measurement_matrix.t() + this->measurement_noise;
  arma::mat kalman_gain = predicted_covariance*this->measurement_matrix.t()*arma::inv(innovation_covariance);
  arma::colvec updated_mean = predicted_mean + kalman_gain*innovation;
  arma::mat eye_mat;
  eye_mat.eye(arma::size(predicted_covariance));
  arma::mat updated_covariance = (eye_mat - kalman_gain*this->measurement_matrix)*predicted_covariance;
  current_state->set_current_posterior_statistics(updated_mean, updated_covariance);
  current_state->log_likelihood = current_state->log_likelihood + dmvnorm(current_measurement,predicted_measurement,innovation_covariance);
}

void ExactKalmanUpdater::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    this->measurement_matrix = this->measurement_matrix_function(conditioned_on_parameters_in);
    this->measurement_noise = this->measurement_noise_function(conditioned_on_parameters_in);
  }
}
