#include "unscented_kalman_predictor.h"
#include "utils.h"
#include "kalman_filter_output.h"

UnscentedKalmanPredictor::UnscentedKalmanPredictor()
  :KalmanPredictor()
{
  this->w0 = 1.0/3.0;
}

UnscentedKalmanPredictor::UnscentedKalmanPredictor(SimulateTransitionKernelPtr transition_kernel_in,
                                                   const arma::mat &process_noise_in,
                                                   double w0_in)
{
  this->transition_kernel = transition_kernel_in;
  this->process_noise = process_noise_in;
  //this->set_using_time = false;
  this->set_using_parameters = false;
  this->w0 = w0_in;
}

UnscentedKalmanPredictor::UnscentedKalmanPredictor(SimulateTransitionKernelFromParametersPtr transition_kernel_function_in,
                                                   GetMatrixPtr process_noise_function_in,
                                                   double w0_in)
{
  this->transition_kernel_function = transition_kernel_function_in;
  this->process_noise_function = process_noise_function_in;
  //this->set_using_time = true;
  this->set_using_parameters = false;
  this->w0 = w0_in;
}

//GetProcessMatrixFromParametersPtr transition_matrix_parameters_function;
//GetProcessMatrixFromParametersPtr process_noise_parameters_function;

/*
UnscentedKalmanPredictor::UnscentedKalmanPredictor(SimulateTransitionKernelFromTimeParametersPtr transition_kernel_time_parameters_function_in,
                                                   GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function_in,
                                                   double w0_in)
{
  this->transition_kernel_time_parameters_function = transition_kernel_time_parameters_function_in;
  this->process_noise_time_parameters_function = process_noise_time_parameters_function_in;
  this->set_using_time = true;
  this->set_using_parameters = true;
  this->w0 = w0_in;
}
*/

UnscentedKalmanPredictor::~UnscentedKalmanPredictor()
{
  
}

UnscentedKalmanPredictor::UnscentedKalmanPredictor(const UnscentedKalmanPredictor &another)
  :KalmanPredictor(another)
{
  this->make_copy(another);
}

void UnscentedKalmanPredictor::operator=(const UnscentedKalmanPredictor &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  KalmanPredictor::operator=(another);
  this->make_copy(another);
}

KalmanPredictor* UnscentedKalmanPredictor::duplicate(void) const
{
  return( new UnscentedKalmanPredictor(*this));
}

void UnscentedKalmanPredictor::make_copy(const UnscentedKalmanPredictor &another)
{
  this->w0 = another.w0;
  this->transition_kernel = another.transition_kernel;
  //this->transition_kernel_time_function = another.transition_kernel_time_function;
  this->transition_kernel_function = another.transition_kernel_function;
  this->process_noise = another.process_noise;
  this->process_noise_function = another.process_noise_function;
  //this->process_noise_time_parameters_function = another.process_noise_time_parameters_function;
}

void UnscentedKalmanPredictor::predict(KalmanFilterOutput* current_state)
{
  //double time_diff = next_time-current_time;
  
  arma::mat sigma_points = get_sigma_points(current_state->posterior_mean_back(),
                                            current_state->posterior_covariance_back(),
                                            this->w0);
  arma::colvec unscented_weights = get_unscented_weights(current_state->posterior_mean_back(),
                                                         this->w0);
  
  arma::mat transformed_sigma_points = arma::mat(arma::size(sigma_points));
  
  arma::colvec predicted_mean = arma::colvec(sigma_points.n_rows,arma::fill::zeros);
  
  arma::mat predicted_covariance;
  if (this->set_using_parameters)
  {
    for (size_t i=0; i<sigma_points.n_cols; ++i)
    {
      transformed_sigma_points.col(i) = this->transition_kernel_function(sigma_points.col(i),
                                                                         this->conditioned_on_parameters);
      predicted_mean = predicted_mean + unscented_weights[i] * transformed_sigma_points.col(i);
      
      predicted_covariance = this->process_noise_function(this->conditioned_on_parameters);
    }
  }
  else
  {
    for (size_t i=0; i<sigma_points.n_cols; ++i)
    {
      transformed_sigma_points.col(i) = this->transition_kernel(sigma_points.col(i));
      predicted_mean = predicted_mean + unscented_weights[i] * transformed_sigma_points.col(i);
      
      predicted_covariance = this->process_noise;
    }
  }
  
  for (size_t i=0; i<sigma_points.n_cols; ++i)
  {
    arma::colvec x_minus_mu = transformed_sigma_points.col(i) - predicted_mean;
    predicted_covariance = predicted_covariance + unscented_weights[i] * x_minus_mu * x_minus_mu.t();
  }
  current_state->set_current_predicted_statistics(predicted_mean, predicted_covariance);
}

void UnscentedKalmanPredictor::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  /*
  if (this->set_using_parameters && !this->set_using_time)
  {
    this->transition_matrix = this->transition_matrix_parameters_function(conditioned_on_parameters_in);
    this->process_noise = this->process_noise_parameters_function(conditioned_on_parameters_in);
  }
  else if (this->set_using_parameters && this->set_using_time)
  {
    this->conditioned_on_parameters = conditioned_on_parameters_in;
  }
  */
  
  if (this->set_using_parameters)
  {
    this->conditioned_on_parameters = conditioned_on_parameters_in;
  }
}
