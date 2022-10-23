#include "exact_kalman_predictor.h"
#include "utils.h"
#include "kalman_filter_output.h"

ExactKalmanPredictor::ExactKalmanPredictor()
  :KalmanPredictor()
{
}

ExactKalmanPredictor::ExactKalmanPredictor(const arma::mat &transition_matrix_in,
                                           const arma::mat &process_noise_in)
{
  this->transition_matrix = transition_matrix_in;
  this->process_noise = process_noise_in;
  this->set_using_time = false;
  this->set_using_parameters = false;
}

ExactKalmanPredictor::ExactKalmanPredictor(GetProcessMatrixFromTimePtr transition_matrix_time_function_in,
                     GetProcessMatrixFromTimePtr process_noise_time_function_in)
{
  this->transition_matrix_time_function = transition_matrix_time_function_in;
  this->process_noise_time_function = process_noise_time_function_in;
  this->set_using_time = true;
  this->set_using_parameters = false;
}

/*
ExactKalmanPredictor::ExactKalmanPredictor(GetProcessMatrixFromParametersPtr transition_matrix_parameters_function_in,
                     GetProcessMatrixFromParametersPtr process_noise_parameters_function_in)
{
  this->transition_matrix_parameters_function = transition_matrix_parameters_function;
  this->process_noise_parameters_function = process_noise_parameters_function;
  this->set_using_time = false;
  this->set_using_parameters = true;
}
 */

ExactKalmanPredictor::ExactKalmanPredictor(GetProcessMatrixFromTimeParametersPtr transition_matrix_time_parameters_function_in,
                     GetProcessMatrixFromTimeParametersPtr process_noise_time_parameters_function_in)
{
  this->transition_matrix_time_parameters_function = transition_matrix_time_parameters_function_in;
  this->process_noise_time_parameters_function = process_noise_time_parameters_function_in;
  this->set_using_time = true;
  this->set_using_parameters = true;
}

ExactKalmanPredictor::~ExactKalmanPredictor()
{
  
}

ExactKalmanPredictor::ExactKalmanPredictor(const ExactKalmanPredictor &another)
  :KalmanPredictor(another)
{
  this->make_copy(another);
}

void ExactKalmanPredictor::operator=(const ExactKalmanPredictor &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  KalmanPredictor::operator=(another);
  this->make_copy(another);
}

KalmanPredictor* ExactKalmanPredictor::duplicate(void) const
{
  return( new ExactKalmanPredictor(*this));
}

void ExactKalmanPredictor::make_copy(const ExactKalmanPredictor &another)
{
  this->transition_matrix = another.transition_matrix;
  this->process_noise = another.process_noise;
  
  this->transition_matrix_time_function = another.transition_matrix_time_function;
  this->process_noise_time_function = another.process_noise_time_function;
  
  //this->transition_matrix_parameters_function = another.transition_matrix_parameters_function;
  //this->process_noise_parameters_function = another.process_noise_parameters_function;
  
  this->transition_matrix_time_parameters_function = another.transition_matrix_time_parameters_function;
  this->process_noise_time_parameters_function = another.process_noise_time_parameters_function;
}

void ExactKalmanPredictor::predict(KalmanFilterOutput* current_state,
                                   double current_time,
                                   double next_time)
{
  double time_diff = next_time-current_time;
  
  if ((!this->set_using_parameters) && this->set_using_time)
  {
    this->transition_matrix = this->transition_matrix_time_function(time_diff);
    this->process_noise = this->process_noise_time_function(time_diff);
  }
  else if (this->set_using_parameters)
  {
    this->transition_matrix = this->transition_matrix_time_parameters_function(time_diff,
                                                                               this->conditioned_on_parameters);
    this->process_noise = this->process_noise_time_parameters_function(time_diff,
                                                                       this->conditioned_on_parameters);
  }
  
  arma::colvec predicted_mean = this->transition_matrix*current_state->posterior_mean_back();
  arma::mat predicted_covariance = this->transition_matrix*current_state->posterior_covariance_back()*this->transition_matrix.t() + this->process_noise;
  current_state->set_current_predicted_statistics(predicted_mean, predicted_covariance);
}

void ExactKalmanPredictor::set_parameters(const Parameters &conditioned_on_parameters_in)
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
