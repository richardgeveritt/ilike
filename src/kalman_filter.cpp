#include "kalman_filter.h"
#include "caller_tempering.h"
#include "kalman_filter_output.h"
#include "kalman_updater.h"
#include "kalman_predictor.h"

namespace ilike
{
KalmanFilter::KalmanFilter()
:LikelihoodEstimator()
{
}

KalmanFilter::KalmanFilter(Data* data_in,
                           size_t lag_in,
                           //const std::string &state_name_in,
                           const arma::colvec &prior_mean_in,
                           const arma::mat &prior_covariance_in,
                           const std::string &index_name_in,
                           const std::string &time_name_in,
                           const std::string &time_diff_name_in,
                           //const std::vector<std::string> &measurements_names_in,
                           size_t first_index_in,
                           size_t last_index_in,
                           size_t predictions_per_update_in,
                           double update_time_step_in,
                           double initial_time_in,
                           bool last_index_is_fixed_in,
                           KalmanPredictor* predictor_in,
                           const std::vector<KalmanUpdater*> &updaters_in,
                           bool smcfixed_flag_in,
                           const std::string &results_name_in)
:LikelihoodEstimator(NULL, NULL, data_in, Parameters(), smcfixed_flag_in)
{
  this->lag = lag_in;
  this->prior_mean = prior_mean_in;
  this->prior_covariance = prior_covariance_in;
  this->state_dimension = prior_mean_in.n_elem;

  this->state_name = "";
  
  for (auto i=updaters_in.begin();
       i!=updaters_in.end();
       ++i)
  {
    if (this->state_name=="")
      this->state_name = (*i)->get_state_variable();
    else
    {
      if (this->state_name != (*i)->get_state_variable())
        Rcpp::stop("KalmanFilter - state variables in all Kalman updaters must be the same.");
    }
  }
  
  this->index_name = index_name_in;
  this->time_name = time_name_in;
  this->time_diff_name = time_diff_name_in;
  this->first_index = first_index_in;
  this->last_index = last_index_in;
  this->predictions_per_update = predictions_per_update_in;
  this->update_time_step = update_time_step_in;
  this->initial_time = initial_time_in;
  this->current_time = initial_time_in;
  this->current_index = this->first_index;
  this->last_index_is_fixed = last_index_is_fixed_in;
  // takes ownership
  this->updaters.clear();
  for (auto i=updaters_in.begin(); i!=updaters_in.end(); ++i)
    this->updaters.emplace_back(*i);
  this->predictor.reset(predictor_in);
  this->results_name = results_name_in;
}

KalmanFilter::~KalmanFilter()
{
  
  this->predictor.reset();
}

//Copy constructor for the KalmanFilter class.
KalmanFilter::KalmanFilter(const KalmanFilter &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void KalmanFilter::operator=(const KalmanFilter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  this->updaters.clear();
  
  this->predictor.reset();
  
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

//LikelihoodEstimator* KalmanFilter::duplicate(void)const
//{
//  return( new KalmanFilter(*this));
//}

void KalmanFilter::make_copy(const KalmanFilter &another)
{
  this->prior_mean = another.prior_mean;
  // B31: the covariance (the copy took the prior mean as its covariance)
  this->prior_covariance = another.prior_covariance;
  this->state_name = another.state_name;
  this->state_dimension = another.state_dimension;
  this->index_name = another.index_name;
  this->time_name = another.time_name;
  this->time_diff_name = another.time_diff_name;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->initial_time = another.initial_time;
  this->current_time = another.current_time;
  this->steady_state_tolerance = another.steady_state_tolerance;
  this->in_steady_state = false;
  this->current_index = another.current_index;
  this->last_index_is_fixed = another.last_index_is_fixed;
  this->lag = another.lag;
  this->results_name = another.results_name;
  this->follows_caller_target = another.follows_caller_target;
  this->followed_variable = another.followed_variable;
  this->followed_anneal = another.followed_anneal;
  this->followed_steps = another.followed_steps;
  
  // its own matrices (the filter writes its time variables in place)
  this->schedule_parameters = another.schedule_parameters.deep_copy();
  
  for (auto i=another.updaters.begin();
       i!=another.updaters.end();
       ++i)
  {
    if (*i!=NULL)
      this->updaters.emplace_back((*i)->duplicate());
    else
      this->updaters.emplace_back(nullptr);
  }
  
  if (another.predictor!=NULL)
    this->predictor.reset(another.predictor->duplicate());
  else
    this->predictor.reset();
}

// double KalmanFilter::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimator* KalmanFilter::duplicate() const
{
  return( new KalmanFilter(*this));
}

void KalmanFilter::setup()
{
  
}

void KalmanFilter::setup(const Parameters &parameters)
{
  
}

KalmanFilterOutput* KalmanFilter::run()
{
  KalmanFilterOutput* current_state = this->kalman_filter_initialise();
  this->evaluate(current_state);
  return current_state;
}

LikelihoodEstimatorOutput* KalmanFilter::initialise()
{
  return this->kalman_filter_initialise();
}

KalmanFilterOutput* KalmanFilter::kalman_filter_initialise()
{
  this->current_index = this->first_index;
  return new KalmanFilterOutput(this,
                                this->lag,
                                this->results_name);
}

void KalmanFilter::evaluate(KalmanFilterOutput* current_state)
{
  if (!this->last_index_is_fixed)
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  this->current_index = this->first_index;
  this->current_time = this->initial_time;
  this->run_steps(current_state, Parameters());
}

void KalmanFilter::subsample_evaluate(KalmanFilterOutput* current_state)
{
  Rcpp::stop("KalmanFilter::subsample_evaluate - not written.");
  
}

void KalmanFilter::subsample_evaluate(KalmanFilterOutput* current_state,
                                      const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("KalmanFilter::subsample_evaluate - not written.");
}

KalmanFilterOutput* KalmanFilter::run(const Parameters &conditioned_on_parameters)
{
  KalmanFilterOutput* current_state = this->kalman_filter_initialise(conditioned_on_parameters);
  this->evaluate(current_state, conditioned_on_parameters);
  return current_state;
}

LikelihoodEstimatorOutput* KalmanFilter::initialise(const Parameters &parameters)
{
  return new KalmanFilterOutput(this,
                                this->lag,
                                this->results_name);
}

KalmanFilterOutput* KalmanFilter::kalman_filter_initialise(const Parameters &parameters)
{
  return new KalmanFilterOutput(this,
                                this->lag,
                                this->results_name);
}

void KalmanFilter::evaluate(KalmanFilterOutput* current_state,
                            const Parameters &conditioned_on_parameters)
{
  if (this->follows_caller_target)
  {
    this->follow_evaluate(current_state, conditioned_on_parameters);
    return;
  }

  // use conditioned_on_parameters to set the next index to stop on
  if (!this->last_index_is_fixed)
  {
    this->first_index = this->last_index+1;
    this->current_index = this->first_index;
    this->last_index = size_t(conditioned_on_parameters[this->index_name][0]);
    if (this->first_index>this->last_index)
      Rcpp::stop("KalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
    this->current_time = this->current_time + this->update_time_step; // (continues from the last run's last step)
  }
  else
  {
    this->current_index = this->first_index;
    this->current_time = this->initial_time;
  }
  this->run_steps(current_state, conditioned_on_parameters);
}

// (C273) the two evaluate functions used to carry two copies of this loop,
// setting the matrices once per run (with the first index, and a time
// difference of 0 in one and the prediction step in the other)
void KalmanFilter::run_steps(KalmanFilterOutput* current_state,
                             const Parameters &conditioned_on_parameters)
{
  if (current_state->predicted_size()==0)
    current_state->start_time = std::chrono::high_resolution_clock::now();
  Parameters schedule, all_parameters;
  this->step_parameters_for_run(conditioned_on_parameters, schedule, all_parameters);
  this->in_steady_state = false; // (C276: a run starts from its prior)
  this->previous_predicted_covariance.reset();
  while (1)
  {
    this->kf_step(current_state, all_parameters, schedule, this->current_index, this->current_time,
                  current_state->predicted_size()==0);
    current_state->set_time();
    if (current_state->results_name!="")
      current_state->write(results_name);
    current_state->start_time = std::chrono::high_resolution_clock::now();
    if (this->current_index >= this->last_index)
    {
      current_state->terminate();
      break;
    }
    current_state->increment_kf_iteration();
    this->current_index = this->current_index + 1;
    this->current_time = this->current_time + this->update_time_step;
  }
}

void KalmanFilter::step_parameters_for_run(const Parameters &conditioned_on_parameters,
                                           Parameters &schedule, Parameters &all_parameters) const
{
  schedule[this->index_name] = 0.0;
  schedule[this->time_name] = 0.0;
  schedule[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  all_parameters.merge_with_fixed(conditioned_on_parameters);
  all_parameters.merge_with_fixed(schedule);
}

void KalmanFilter::kf_step(KalmanFilterOutput* current_state,
                           const Parameters &all_parameters,
                           Parameters &schedule, size_t index, double time, bool initial)
{
  schedule[this->index_name] = double(index); // (seen through all_parameters)
  schedule[this->time_name] = time;

  // the matrices for this step (the k predictions all see the step's values,
  // as the particle filter's k moves do)
  if (!initial)
    this->predictor->set_parameters(all_parameters);
  for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
    (*i)->set_parameters(all_parameters);

  // (C276) the steady-state option
  bool steady = (this->steady_state_tolerance > 0.0) && !this->follows_caller_target && !initial;
  bool unchanged = false;
  if (steady)
  {
    unchanged = !this->predictor->matrices_changed();
    for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
      if ((*i)->matrices_changed())
        unchanged = false;
    if (!unchanged)
      this->in_steady_state = false;
  }
  if (steady && this->in_steady_state)
  {
    // the mean alone, with the frozen covariances and gains
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
      this->predictor->predict_mean(current_state);
    current_state->set_current_predicted_covariance(this->frozen_predicted_covariance);
    current_state->add_predicted_statistics();
    current_state->set_current_posterior_to_be_current_predicted();
    for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
    {
      arma::colvec current_measurement = arma::vectorise(this->data->row(index)[(*i)->get_measurement_variable()]);
      (*i)->update_mean(current_state, current_measurement);
    }
    current_state->set_current_posterior_covariance(this->frozen_posterior_covariance);
    current_state->set_llhd(current_state->log_likelihood);
    current_state->add_posterior_statistics();
    current_state->add_schedule_parameters(schedule.deep_copy());
    return;
  }

  if (initial)
  {
    current_state->set_current_predicted_statistics(this->prior_mean, this->prior_covariance);
  }
  else
  {
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
      this->predictor->predict(current_state);
  }
  current_state->add_predicted_statistics();

  // (C276) converged: the matrices are the same as at the last step, and the
  // predicted covariance has changed by at most the tolerance; the gains are
  // frozen in this step's updates
  bool freeze = false;
  if (steady)
  {
    const arma::mat &P = current_state->get_current_predicted_covariance();
    if (unchanged && (this->previous_predicted_covariance.n_elem == P.n_elem) && (P.n_elem > 0))
    {
      double change = arma::abs(P - this->previous_predicted_covariance).max();
      double size = arma::abs(P).max();
      freeze = (change <= this->steady_state_tolerance*size);
    }
    this->previous_predicted_covariance = P;
  }

  current_state->set_current_posterior_to_be_current_predicted(); // (C271)
  for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
  {
    if (freeze)
      (*i)->freeze_gain(current_state->get_current_posterior_covariance());
    arma::colvec current_measurement = arma::vectorise(this->data->row(index)[(*i)->get_measurement_variable()]);
    (*i)->update(current_state, current_measurement);
  }
  if (freeze)
  {
    this->frozen_predicted_covariance = current_state->get_current_predicted_covariance();
    this->frozen_posterior_covariance = current_state->get_current_posterior_covariance();
    this->in_steady_state = true;
  }
  current_state->set_llhd(current_state->log_likelihood);
  current_state->add_posterior_statistics();
  current_state->add_schedule_parameters(schedule.deep_copy());
}

void KalmanFilter::set_steady_state_tolerance(double tolerance)
{
  if (!(tolerance >= 0.0))
    Rcpp::stop("kf: steady_state must be a tolerance of 0 or more.");
  this->steady_state_tolerance = tolerance;
  if (tolerance > 0.0)
  {
    if (!this->steady_state_possible())
      Rcpp::stop("kf: steady_state is available for linear-Gaussian transition and data models only.");
    this->predictor->track_matrix_changes();
    for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
      (*i)->track_matrix_changes();
  }
}

bool KalmanFilter::steady_state_possible() const
{
  if ((this->predictor == NULL) || !this->predictor->supports_steady_state())
    return false;
  for (auto i=this->updaters.begin(); i!=this->updaters.end(); ++i)
    if (!(*i)->supports_steady_state())
      return false;
  return true;
}

bool KalmanFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}

void KalmanFilter::specific_change_data(Data* new_data)
{
  
}

// --- SMC^2 with a Kalman filter (dev/plan.md 3.21 (c2''')) --------------------

void KalmanFilter::set_follows_caller_target(const std::string &fraction_variable, bool anneal)
{
  this->follows_caller_target = true;
  this->followed_variable = fraction_variable;
  this->followed_anneal = anneal;
  this->smcfixed_flag = false; // (the value changes along the caller's sequence)
}

// one predict-and-update step: the output's next observation (the first one:
// the prior, then the update, as a run's initial step does), with the step's
// index, time and time difference (C273: kf_step, shared with evaluate)
void KalmanFilter::follow_step(KalmanFilterOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  size_t step = current_state->followed_steps; // (observations already in this output)
  Parameters schedule, all_parameters;
  this->step_parameters_for_run(conditioned_on_parameters, schedule, all_parameters);
  this->kf_step(current_state, all_parameters, schedule, this->first_index + step,
                this->initial_time + double(step)*this->update_time_step, step == 0);
  current_state->set_time(); // (the output's times are written beside its other records)
  current_state->followed_steps = step + 1;
}

// keep the next observation's update: the output's log-likelihood starts
// from the committed value (a trial's value is not a base for the update)
void KalmanFilter::follow_commit(KalmanFilterOutput* current_state,
                                 const Parameters &conditioned_on_parameters)
{
  current_state->log_likelihood = current_state->committed_log_likelihood;
  this->follow_step(current_state, conditioned_on_parameters);
  current_state->committed_log_likelihood = current_state->log_likelihood;
  if (current_state->followed_steps > this->followed_steps)
    this->followed_steps = current_state->followed_steps;
}

// the value at the caller's target (k, r): the committed observations, plus
// the next one's exact log predictive density (times r with anneal), from a
// step on a copy of the output; no random numbers
double KalmanFilter::follow_trial_value(KalmanFilterOutput* current_state,
                                        const Parameters &conditioned_on_parameters,
                                        size_t k, double r)
{
  double value = current_state->committed_log_likelihood;
  if (value == -arma::datum::inf)
    return value;
  switch (caller_trial(current_state->followed_steps, k, r))
  {
    case CallerTrial::committed:
      return value;
    case CallerTrial::next:
    {
      KalmanFilterOutput trial(*current_state);
      trial.log_likelihood = value;
      this->follow_step(&trial, conditioned_on_parameters);
      return value + ((r > 0.0) ? r : 1.0)*(trial.log_likelihood - value);
    }
    default:
      Rcpp::stop(caller_trial_too_far_message());
  }
  return value;
}

// a new output: run through the observations committed so far (the same
// steps as an output that has followed the caller from its start), then a
// trial at the caller's current target
void KalmanFilter::follow_from_start(KalmanFilterOutput* current_state,
                                     const Parameters &conditioned_on_parameters)
{
  size_t k;
  double r;
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  size_t committed = std::min(this->followed_steps, k);
  current_state->followed_steps = 0;
  current_state->committed_log_likelihood = 0.0;
  current_state->log_likelihood = 0.0;
  while ((current_state->followed_steps < committed) && (current_state->committed_log_likelihood != -arma::datum::inf))
    this->follow_commit(current_state, conditioned_on_parameters);
  current_state->caught_up = true;
  current_state->log_likelihood = this->follow_trial_value(current_state, conditioned_on_parameters, k, r);
}

// an evaluation: a new output is run from the start; afterwards a trial at
// the caller's target
void KalmanFilter::follow_evaluate(KalmanFilterOutput* current_state,
                                   const Parameters &conditioned_on_parameters)
{
  if (!current_state->caught_up)
  {
    this->follow_from_start(current_state, conditioned_on_parameters);
    return;
  }
  size_t k;
  double r;
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  current_state->log_likelihood = this->follow_trial_value(current_state, conditioned_on_parameters, k, r);
}

// the caller's move stage: commit the next observation once it is whole
// (nothing random to do: the step is deterministic); the value stays the one
// at the caller's current target
void KalmanFilter::advance_at_caller_target_kf(KalmanFilterOutput* current_state,
                                               const Parameters &conditioned_on_parameters)
{
  if (!this->follows_caller_target || !current_state->caught_up)
    return;
  size_t k;
  double r;
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  if ((current_state->committed_log_likelihood != -arma::datum::inf) && (k == current_state->followed_steps + 1) && (r == 0.0))
    this->follow_commit(current_state, conditioned_on_parameters);
  current_state->log_likelihood = this->follow_trial_value(current_state, conditioned_on_parameters, k, r);
}

}
