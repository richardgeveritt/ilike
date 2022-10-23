#include "kalman_filter.h"
#include "kalman_filter_output.h"
#include "kalman_updater.h"
#include "kalman_predictor.h"

KalmanFilter::KalmanFilter()
  :LikelihoodEstimator()
{
}

KalmanFilter::KalmanFilter(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           Data* data_in,
                           size_t lag_in,
                           EvaluateLogLikelihoodPtr llhd_in,
                           double current_time_in,
                           bool smcfixed_flag_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->current_time = current_time_in;
  this->smcfixed_flag = smcfixed_flag_in;
  this->lag = lag_in;
  //this->output = new KalmanFilterOutput();
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
  this->smcfixed_flag = another.smcfixed_flag;
  this->prior_mean = another.prior_mean;
  this->prior_covariance = another.prior_mean;
  this->index_name = another.index_name;
  this->measurements_names = another.measurements_names;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->current_time = another.current_time;
  this->current_index = another.current_index;
  this->last_index_is_fixed = another.last_index_is_fixed;
  this->lag = another.lag;
  
  if (another.updater!=NULL)
    this->updater = another.updater->duplicate();
  
  if (another.predictor!=NULL)
    this->predictor = another.predictor->duplicate();
}

// double KalmanFilter::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }


KalmanFilterOutput* KalmanFilter::run()
{
  KalmanFilterOutput* current_state = this->kalman_filter_initialise();
  this->evaluate(current_state);
  return current_state;
}

LikelihoodEstimatorOutput* KalmanFilter::initialise()
{
  this->first_index = 0;
  this->current_index = this->first_index;
  return new KalmanFilterOutput(this,
                                this->lag);
}

KalmanFilterOutput* KalmanFilter::kalman_filter_initialise()
{
  this->first_index = 0;
  this->current_index = this->first_index;
  return new KalmanFilterOutput(this,
                                this->lag);
}

void KalmanFilter::evaluate(KalmanFilterOutput* current_state)
{
  
  if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  }
  
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
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
    arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
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
    
    arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
  }
}

void KalmanFilter::subsample_evaluate(KalmanFilterOutput* current_state)
{
  throw std::runtime_error("KalmanFilter::subsample_evaluate - not written.");
  
  /*
  if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  }
  
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
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
  throw std::runtime_error("KalmanFilter::subsample_evaluate - not written.");
  /*
  if ( (this->updater->set_using_parameters) || (this->updater->set_using_parameters) )
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter fix updater and/or predictor.");
  }
  
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("KalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
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
                                this->lag);
}

KalmanFilterOutput* KalmanFilter::kalman_filter_initialise(const Parameters &parameters)
{
  return new KalmanFilterOutput(this,
                                this->lag);
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
      throw std::runtime_error("KalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
  }
  
  // Set predictor and updater with parameters.
  this->predictor->set_parameters(conditioned_on_parameters);
  this->updater->set_parameters(conditioned_on_parameters);
  
  if (current_state->predicted_size()==0)
  {
    // Initial step.
    current_state->set_current_predicted_statistics(this->prior_mean,
                                                    this->prior_covariance);
    current_state->add_predicted_statistics();
    
    // For a particle filter, we instead need to use
    // Data current_measurement = this->data->get_using_time_index(this->measurements_names);
    // Returns the data for a time slice, which will include a bunch of variables.
    arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
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
    
    arma::colvec current_measurement = (*this->data)[this->measurements_names].col(this->current_index);
    
    this->updater->update(current_state,
                          current_measurement);
    current_state->add_posterior_statistics();
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
