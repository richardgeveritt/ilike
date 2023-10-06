#include "unscented_kalman_updater.h"
#include "utils.h"
#include "kalman_filter_output.h"

#include <math.h>

UnscentedKalmanUpdater::UnscentedKalmanUpdater()
  :KalmanUpdater()
{
  this->w0 = 1.0/3.0;
}

UnscentedKalmanUpdater::UnscentedKalmanUpdater(GetMatrixSimulateMeasurementKernelPtr measurement_kernel_function_in,
                                               GetMeasurementMatrixPtr measurement_noise_function_in,
                                               double w0_in)
{
  this->measurement_kernel_function = measurement_kernel_function_in;
  this->measurement_noise_function = measurement_noise_function_in;
  this->w0 = w0_in;
}

UnscentedKalmanUpdater::UnscentedKalmanUpdater(MatrixSimulateMeasurementKernelPtr measurement_kernel_in,
                                               const arma::mat &measurement_noise_in,
                                               double w0_in)
{
  this->measurement_kernel = measurement_kernel_in;
  this->measurement_noise = measurement_noise_in;
  this->w0 = w0_in;
}

UnscentedKalmanUpdater::~UnscentedKalmanUpdater()
{
  
}

UnscentedKalmanUpdater::UnscentedKalmanUpdater(const UnscentedKalmanUpdater &another)
  :KalmanUpdater(another)
{
  this->make_copy(another);
}

void UnscentedKalmanUpdater::operator=(const UnscentedKalmanUpdater &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  KalmanUpdater::operator=(another);
  this->make_copy(another);
}

KalmanUpdater* UnscentedKalmanUpdater::duplicate() const
{
  return( new UnscentedKalmanUpdater(*this));
}

void UnscentedKalmanUpdater::make_copy(const UnscentedKalmanUpdater &another)
{
  this->w0 = another.w0;
  this->measurement_kernel = another.measurement_kernel;
  this->measurement_noise = another.measurement_noise;
  this->measurement_kernel_function = another.measurement_kernel_function;
  this->measurement_noise_function = another.measurement_noise_function;
  this->conditioned_on_parameters = another.conditioned_on_parameters;
}

void UnscentedKalmanUpdater::update(KalmanFilterOutput* current_state,
                                    const arma::colvec &current_measurement)
{
  arma::colvec predicted_mean = current_state->predicted_mean_back();
  arma::colvec predicted_covariance = current_state->predicted_covariance_back();
  arma::mat sigma_points = get_sigma_points(predicted_mean,
                                            predicted_covariance,
                                            this->w0);
  arma::colvec unscented_weights = get_unscented_weights(predicted_mean,
                                                         this->w0);
  
  arma::mat transformed_sigma_points = arma::mat(current_measurement.n_rows,
                                                 sigma_points.n_cols);
  
  arma::colvec predicted_measurement = arma::colvec(current_measurement.n_rows,arma::fill::zeros);
  
  if (this->set_using_parameters)
  {
    for (size_t i=0; i<sigma_points.n_cols; ++i)
    {
      transformed_sigma_points.col(i) = this->measurement_kernel_function(sigma_points.col(i),
                                                                          this->conditioned_on_parameters);
      predicted_measurement = predicted_measurement + unscented_weights[i] * transformed_sigma_points.col(i);
    }
  }
  else
  {
    for (size_t i=0; i<sigma_points.n_cols; ++i)
    {
      transformed_sigma_points.col(i) = this->measurement_kernel(sigma_points.col(i));
      predicted_measurement = predicted_measurement + unscented_weights[i] * transformed_sigma_points.col(i);
    }
  }
  
  arma::mat innovation_covariance = this->measurement_noise;
  arma::mat cross_covariance = arma::colvec(sigma_points.n_rows,
                                            current_measurement.n_rows,
                                            arma::fill::zeros);
  for (size_t i=0; i<sigma_points.n_cols; ++i)
  {
    arma::colvec x_minus_mu = transformed_sigma_points.col(i) - predicted_measurement;
    innovation_covariance = innovation_covariance + unscented_weights[i] * x_minus_mu * x_minus_mu.t();
    cross_covariance = cross_covariance + unscented_weights[i] * (sigma_points.col(i) - predicted_mean) * x_minus_mu.t();
  }
  
  arma::colvec innovation = current_measurement - predicted_measurement;
  arma::mat kalman_gain = cross_covariance*arma::inv(innovation_covariance);
  arma::colvec updated_mean = predicted_mean + kalman_gain*innovation;
  arma::mat updated_covariance = predicted_covariance - kalman_gain*innovation_covariance*kalman_gain.t();
  current_state->set_current_posterior_statistics(updated_mean, updated_covariance);
  current_state->log_likelihood = current_state->log_likelihood + dmvnorm(current_measurement,predicted_measurement,innovation_covariance);
}

void UnscentedKalmanUpdater::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    this->conditioned_on_parameters = conditioned_on_parameters_in;
    this->measurement_noise = this->measurement_noise_function(conditioned_on_parameters_in);
  }
}
