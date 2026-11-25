#include "exact_kalman_predictor.h"
#include "utils.h"
#include "kalman_filter_output.h"

namespace ilike
{

// (C276) exactly the same matrix (sizes and every value)
static bool same_matrix(const arma::mat &a, const arma::mat &b)
{
  return (a.n_rows == b.n_rows) && (a.n_cols == b.n_cols) &&
    std::equal(a.memptr(), a.memptr() + a.n_elem, b.memptr());
}
ExactKalmanPredictor::ExactKalmanPredictor()
:KalmanPredictor()
{
}

ExactKalmanPredictor::ExactKalmanPredictor(const arma::mat &transition_matrix_in,
                                           const arma::mat &process_noise_in)
{
  this->transition_matrix = transition_matrix_in;
  this->process_noise = process_noise_in;
  this->set_using_parameters = false;
}

ExactKalmanPredictor::ExactKalmanPredictor(GetMatrixPtr transition_matrix_function_in,
                                           GetMatrixPtr process_noise_function_in)
{
  this->transition_matrix_function = transition_matrix_function_in;
  this->process_noise_function = process_noise_function_in;
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

KalmanPredictor* ExactKalmanPredictor::duplicate() const
{
  return( new ExactKalmanPredictor(*this));
}

void ExactKalmanPredictor::make_copy(const ExactKalmanPredictor &another)
{
  this->transition_matrix = another.transition_matrix;
  this->process_noise = another.process_noise;
  
  this->transition_matrix_function = another.transition_matrix_function;
  this->process_noise_function = another.process_noise_function;
  this->tracks_matrix_changes = another.tracks_matrix_changes; // (C276)
  this->changed_at_last_set = another.changed_at_last_set;

}

// (C271) one prediction from the current prediction, which the filter sets
// to the posterior before the first: predictions_per_update > 1 chains them
// (each used to start again from the stored posterior). The matrices are
// those set_parameters evaluated.
void ExactKalmanPredictor::predict(KalmanFilterOutput* current_state)
{
  // (C274) into the workspace members (memory reused between steps)
  this->new_mean = this->transition_matrix*current_state->get_current_predicted_mean();
  this->new_covariance = this->transition_matrix*current_state->get_current_predicted_covariance()*this->transition_matrix.t() + this->process_noise;
  current_state->set_current_predicted_statistics(this->new_mean, this->new_covariance);
}

// (C271, C273) the matrix functions are evaluated here, once for each step
// (the filter sets the step's index, time and time difference), instead of
// at every one of the step's predictions
void ExactKalmanPredictor::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    if (this->tracks_matrix_changes)
    {
      arma::mat A = this->transition_matrix_function(conditioned_on_parameters_in);
      arma::mat Q = this->process_noise_function(conditioned_on_parameters_in);
      this->changed_at_last_set = !same_matrix(A, this->transition_matrix) || !same_matrix(Q, this->process_noise);
      this->transition_matrix = std::move(A);
      this->process_noise = std::move(Q);
    }
    else
    {
      this->transition_matrix = this->transition_matrix_function(conditioned_on_parameters_in);
      this->process_noise = this->process_noise_function(conditioned_on_parameters_in);
    }
  }
}

// (C276) in steady state: the mean alone (the filter keeps the covariance)
void ExactKalmanPredictor::predict_mean(KalmanFilterOutput* current_state)
{
  this->new_mean = this->transition_matrix*current_state->get_current_predicted_mean();
  current_state->set_current_predicted_mean(this->new_mean);
}
}
