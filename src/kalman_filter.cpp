#include "kalman_filter.h"
#include "kalman_filter_output.h"
#include "kalman_updater.h"
#include "kalman_predictor.h"

KalmanFilter::KalmanFilter()
  :LikelihoodEstimator()
{
}

KalmanFilter::KalmanFilter(Data* data_in,
                           size_t lag_in,
                           const std::string &state_name_in,
                           const arma::colvec &prior_mean_in,
                           const arma::mat &prior_covariance_in,
                           const std::string &index_name_in,
                           const std::string &time_name_in,
                           const std::string &time_diff_name_in,
                           const std::vector<std::string> &measurements_names_in,
                           size_t first_index_in,
                           size_t last_index_in,
                           size_t predictions_per_update_in,
                           double update_time_step_in,
                           double initial_time_in,
                           bool last_index_is_fixed_in,
                           KalmanPredictor* predictor_in,
                           KalmanUpdater* updater_in,
                           bool smcfixed_flag_in,
                           const std::string &results_name_in)
:LikelihoodEstimator(NULL, NULL, data_in, Parameters(), smcfixed_flag_in)
{
  this->lag = lag_in;
  this->prior_mean = prior_mean_in;
  this->prior_covariance = prior_covariance_in;
  this->state_dimension = prior_mean_in.n_elem;
  
  this->state_name = state_name_in;
  this->index_name = index_name_in;
  this->time_name = time_name_in;
  this->time_diff_name = time_diff_name_in;
  this->measurements_names = measurements_names_in;
  this->first_index = first_index_in;
  this->last_index = last_index_in;
  this->predictions_per_update = predictions_per_update_in;
  this->update_time_step = update_time_step_in;
  this->current_time = initial_time_in;
  this->current_index = this->first_index;
  this->last_index_is_fixed = last_index_is_fixed_in;
  this->updater = updater_in;
  this->predictor = predictor_in;
  this->results_name = results_name_in;
}

KalmanFilter::~KalmanFilter()
{
  if (this->updater!=NULL)
    delete this->updater;
  
  if (this->predictor!=NULL)
    delete this->predictor;
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
  
  if (this->updater!=NULL)
    delete this->updater;
  
  if (this->predictor!=NULL)
    delete this->predictor;

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
  this->prior_covariance = another.prior_mean;
  this->state_name = another.state_name;
  this->state_dimension = another.state_dimension;
  this->index_name = another.index_name;
  this->time_name = another.time_name;
  this->time_diff_name = another.time_diff_name;
  this->measurements_names = another.measurements_names;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->current_time = another.current_time;
  this->current_index = another.current_index;
  this->last_index_is_fixed = another.last_index_is_fixed;
  this->lag = another.lag;
  this->results_name = another.results_name;
  
  this->schedule_parameters = another.schedule_parameters;
  
  if (another.updater!=NULL)
    this->updater = another.updater->duplicate();
  else
    this->updater = NULL;
  
  if (another.predictor!=NULL)
    this->predictor = another.predictor->duplicate();
  else
    this->predictor = NULL;
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
  //this->first_index = 0;
  this->current_index = this->first_index;
  return new KalmanFilterOutput(this,
                                this->lag,
                                this->results_name);
}

void KalmanFilter::evaluate(KalmanFilterOutput* current_state)
{
  
  //if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  //{
  //  Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  //}
  
  if (!this->last_index_is_fixed)
  {
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->schedule_parameters = Parameters(this->time_name,this->current_time);
  this->schedule_parameters[this->index_name] = this->current_index;
  this->schedule_parameters[this->time_diff_name] = 0.0;
  
  Parameters all_parameters;
  all_parameters.merge_with_fixed(this->schedule_parameters);
  
  // Set predictor and updater with parameters.
  this->predictor->set_parameters(all_parameters);
  this->updater->set_parameters(all_parameters);
  
  if (current_state->predicted_size()==0)
  {
    // Initial step.
    current_state->set_current_predicted_statistics(this->prior_mean,
                                                    this->prior_covariance);
    current_state->add_predicted_statistics();
    
    // For a particle filter, we instead need to use
    // Data current_measurement = this->data->get_using_time_index(this->measurements_names);
    // Returns the data for a time slice, which will include a bunch of variables.
    arma::colvec current_measurement = this->data->row(this->current_index)[this->measurements_names];
    //arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
    // Update at initial step.
    this->updater->update(current_state,
                          current_measurement);
    
    current_state->llhds.push_back(current_state->log_likelihood);
    
    current_state->add_posterior_statistics();
    current_state->add_schedule_parameters(this->schedule_parameters);
    
    current_state->set_time();
    if (current_state->results_name!="")
      current_state->write(results_name);
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    current_state->increment_kf_iteration();
  }
  
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);
  
  while (!this->check_termination())
  {
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
    {
      this->current_time = this->current_time + predict_time_step;
      
      this->schedule_parameters[this->time_name] = this->current_time;
      this->schedule_parameters[this->time_diff_name] = predict_time_step;
      
      //this->predictor->predict(current_state,
      //                         previous_time,
      //                         this->current_time);
      
      this->predictor->predict(current_state);
    }
    current_state->add_predicted_statistics();
    this->current_index = this->current_index + 1;
    
    this->schedule_parameters[this->index_name] = this->current_index;
    
    arma::colvec current_measurement = this->data->row(this->current_index)[this->measurements_names];
    
    //arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    
    current_state->llhds.push_back(current_state->log_likelihood);
    
    current_state->add_posterior_statistics();
    current_state->add_schedule_parameters(this->schedule_parameters);
    
    current_state->set_time();
    if (current_state->results_name!="")
      current_state->write(results_name);
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    current_state->increment_kf_iteration();
  }
}

void KalmanFilter::subsample_evaluate(KalmanFilterOutput* current_state)
{
  Rcpp::stop("KalmanFilter::subsample_evaluate - not written.");
  
  /*
  if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  {
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  }
  
  if (!this->last_index_is_fixed)
  {
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  if (current_state->predicted_size()==0)
  {
    // Initial step.
    current_state->set_current_predicted_statistics(this->prior_mean,
                                                    this->prior_covariance);
    current_state->add_predicted_statistics();
    
    // For a particle filter, we instead need to use
    // Data current_measurement = this->data->get_using_time_index(this->measurements_names);
    // Returns the data for a time slice, which will include a bunch of variables.
    arma::colvec current_measurement = (*this->subsampler->small_data)[this->measurements_names].col(this->current_index);
    
    // Update at initial step.
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
  }
  
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);
  
  while (!this->check_termination())
  {
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
    {
      double previous_time = this->current_time;
      this->current_time = this->current_time + predict_time_step;
      this->predictor->predict(current_state,
                               previous_time,
                               this->current_time);
    }
    current_state->add_predicted_statistics();
    this->current_index = this->current_index + 1;
    
    arma::colvec current_measurement = (*this->subsampler->small_data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
  }
  */
}

void KalmanFilter::subsample_evaluate(KalmanFilterOutput* current_state,
                                      const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("KalmanFilter::subsample_evaluate - not written.");
  /*
  if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  {
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  }
  
  if (!this->last_index_is_fixed)
  {
    Rcpp::stop("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  if (current_state->predicted_size()==0)
  {
    // Initial step.
    current_state->set_current_predicted_statistics(this->prior_mean,
                                                    this->prior_covariance);
    current_state->add_predicted_statistics();
    
    // For a particle filter, we instead need to use
    // Data current_measurement = this->data->get_using_time_index(this->measurements_names);
    // Returns the data for a time slice, which will include a bunch of variables.
    arma::colvec current_measurement = (*this->subsampler->small_data)[this->measurements_names].col(this->current_index);
    
    // Update at initial step.
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
  }
  
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);
  
  while (!this->check_termination())
  {
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
    {
      double previous_time = this->current_time;
      this->current_time = this->current_time + predict_time_step;
      this->predictor->predict(current_state,
                               previous_time,
                               this->current_time);
    }
    current_state->add_predicted_statistics();
    this->current_index = this->current_index + 1;
    
    arma::colvec current_measurement = (*this->subsampler->small_data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
  }
   */
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
  
  // use conditioned_on_parameters to set the next index to stop on
  if (!this->last_index_is_fixed)
  {
    this->first_index = this->last_index+1;
    this->current_index = this->first_index;
    this->last_index = size_t(conditioned_on_parameters[this->index_name][0]);
    if (this->first_index>this->last_index)
      Rcpp::stop("KalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
  }
  
  this->schedule_parameters = Parameters(this->time_name,this->current_time);
  this->schedule_parameters[this->index_name] = this->current_index;
  this->schedule_parameters[this->time_diff_name] = 0.0;
  
  Parameters all_parameters;
  all_parameters.merge_with_fixed(conditioned_on_parameters);
  all_parameters.merge_with_fixed(this->schedule_parameters);
  
  // Set predictor and updater with parameters.
  this->predictor->set_parameters(all_parameters);
  this->updater->set_parameters(all_parameters);
  
  if (current_state->predicted_size()==0)
  {
    // Initial step.
    current_state->set_current_predicted_statistics(this->prior_mean,
                                                    this->prior_covariance);
    current_state->add_predicted_statistics();
    
    // For a particle filter, we instead need to use
    // Data current_measurement = this->data->get_using_time_index(this->measurements_names);
    // Returns the data for a time slice, which will include a bunch of variables.
    arma::colvec current_measurement = this->data->row(this->current_index)[this->measurements_names];
    //arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
    // Update at initial step.
    this->updater->update(current_state,
                          current_measurement);
    
    current_state->llhds.push_back(current_state->log_likelihood);
    
    current_state->add_posterior_statistics();
    current_state->add_schedule_parameters(this->schedule_parameters);
    
    current_state->set_time();
    if (current_state->results_name!="")
      current_state->write(results_name);
    current_state->start_time = std::chrono::high_resolution_clock::now();
  }
  
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);

  while (!this->check_termination())
  {
    current_state->set_current_predicted_to_be_current_posterior();
    for (size_t i=0; i<this->predictions_per_update; ++i)
    {
      //double previous_time = this->current_time;
      this->current_time = this->current_time + predict_time_step;
      
      this->schedule_parameters[this->time_name] = this->current_time;
      this->schedule_parameters[this->time_diff_name] = predict_time_step;
      
      //this->predictor->predict(current_state,
      //                         previous_time,
      //                         this->current_time);
      
      this->predictor->predict(current_state);
    }
    current_state->add_predicted_statistics();
    this->current_index = this->current_index + 1;
    this->schedule_parameters[this->index_name] = this->current_index;
    
    arma::colvec current_measurement = this->data->row(this->current_index)[this->measurements_names];//(*this->data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    
    current_state->llhds.push_back(current_state->log_likelihood);
    
    current_state->add_posterior_statistics();
    current_state->add_schedule_parameters(this->schedule_parameters);
    
    current_state->set_time();
    if (current_state->results_name!="")
      current_state->write(results_name);
    current_state->start_time = std::chrono::high_resolution_clock::now();
  }
}

bool KalmanFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}


//double KalmanFilter::evaluate(const Parameters &parameters)
//{
//  return this->func(parameters,*this->data);
//}

// void KalmanFilter::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
