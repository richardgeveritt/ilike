#include "ensemble_kalman_filter.h"
#include "ensemble_kalman_output.h"
#include "factors.h"
#include "single_point_move_output.h"
#include "proposal_kernel.h"
//#include "ensemble_member.h"
#include "ensemble_factor_variables.h"
#include "ensemble_kalman_filter.h"
#include "ensemble_kalman_worker.h"
#include "ensemble_factors.h"
#include "vector_single_index.h"
#include "single_point_move_output.h"
//#include "ensemble_kalman_updater.h"
//#include "ensemble_kalman_predictor.h"

EnsembleKalmanFilter::EnsembleKalmanFilter()
  :EnsembleKalman()
{
  this->proposal_kernel = NULL;
}

EnsembleKalmanFilter::EnsembleKalmanFilter(RandomNumberGenerator* rng_in,
                                           size_t* seed_in,
                                           Data* data_in,
                                           EvaluateLogLikelihoodPtr llhd_in,
                                           double current_time_in,
                                           bool sequencer_limit_is_fixed_in)
:EnsembleKalman(rng_in, seed_in, data_in, false, true)
{
  this->current_time = current_time_in;
  this->current_index = this->first_index;
  this->proposal_kernel = NULL;
  //this->output = new EnsembleKalmanFilterOutput();
}

EnsembleKalmanFilter::~EnsembleKalmanFilter()
{
  //if (this->proposal_kernel!=NULL)
  //  delete this->proposal_kernel;
}

//Copy constructor for the EnsembleKalmanFilter class.
EnsembleKalmanFilter::EnsembleKalmanFilter(const EnsembleKalmanFilter &another)
  :EnsembleKalman(another)
{
  this->make_copy(another);
}

void EnsembleKalmanFilter::operator=(const EnsembleKalmanFilter &another)
{
  if(this == &another){ //if a==a
    return;
  }

  //if (this->proposal_kernel!=NULL)
  //  delete this->proposal_kernel;
  
  EnsembleKalman::operator=(another);
  this->make_copy(another);
}

//LikelihoodEstimator* EnsembleKalmanFilter::duplicate(void)const
//{
//  return( new EnsembleKalmanFilter(*this));
//}

void EnsembleKalmanFilter::make_copy(const EnsembleKalmanFilter &another)
{
  //if (another.proposal_kernel!=NULL)
  this->proposal_kernel = another.proposal_kernel;
  
  this->index_name = another.index_name;
  //this->measurements_names = another.measurements_names;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->current_time = another.current_time;
  this->current_index = another.current_index;
  this->last_index_is_fixed = another.last_index_is_fixed;
  this->lag = another.lag;
  
  this->smcfixed_flag = another.smcfixed_flag;
  //if (this->output!=NULL)
  //  this->output = another.output->duplicate();
}

// double EnsembleKalmanFilter::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_run()
{
  EnsembleKalmanOutput* current_state = this->ensemble_kalman_initialise();
  this->ensemble_kalman_simulate(current_state);
  this->ensemble_kalman_evaluate(current_state);
  return current_state;
}

EnsembleKalmanOutput* EnsembleKalmanFilter::ensemble_kalman_initialise()
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->factors->set_data(this->current_index);
  return output;
}

void EnsembleKalmanFilter::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state)
{
  if (current_state->all_ensembles.size()==0)
  {
    VectorSingleIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble();
    this->the_worker->move(next_ensemble,
                           current_ensemble);
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state)
{
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("EnsembleKalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->ensemble_kalman_evaluate_smcfixed_part(current_state);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state);
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state)
{
  //this->the_worker->smcfixed_weight(current_state->all_ensembles.back());
  //current_state->initialise_next_step();
}

MoveOutput* EnsembleKalmanFilter::move(RandomNumberGenerator &rng,
                                       Particle &ensemble_member)
{
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    VectorSingleIndex index(this->current_index);
    if (i==0)
    {
      moved_ensemble_member = this->proposal_kernel->move(rng,
                                                          ensemble_member);
    }
    else
    {
      moved_ensemble_member = this->proposal_kernel->move(rng,
                                                          moved_ensemble_member);
    }
  }
  
  moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
    
  MoveOutput* moved_output = new SinglePointMoveOutput(moved_ensemble_member);
  
  return moved_output;
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    std::runtime_error("EnsembleKalmanFilter::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    this->the_worker->shift(&current_state->back(),1.0);
    this->the_worker->unpack(&current_state->back());
    
    /*
    if (current_state->all_particles.size()==1)
    {
      
      this->the_worker->weight(current_state->back());
    }
    else
    {
      if (this->dynamic_proposal_is_evaluated)
      {
        this->the_worker->pf_weight(current_state->back(),
                                    *(current_state->all_particles.end()-2),
                                    this->proposal_kernel);
      }
      else
      {
        this->the_worker->pf_weight(current_state->back(),
                                    *(current_state->all_particles.end()-2),
                                    NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    */
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_simulate(current_state);
    
  }
}

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_run(const Parameters &conditioned_on_parameters)
{
  EnsembleKalmanOutput* current_state = this->ensemble_kalman_initialise(conditioned_on_parameters);
  this->ensemble_kalman_evaluate(current_state, conditioned_on_parameters);
  return current_state;
}

EnsembleKalmanOutput* EnsembleKalmanFilter::ensemble_kalman_initialise(const Parameters &parameters)
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->factors->set_data(this->current_index);
  return output;
}

void EnsembleKalmanFilter::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state,
                                                    const Parameters &conditioned_on_parameters)
{
  if (current_state->all_ensembles.size()==0)
  {
    VectorSingleIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index,
                            conditioned_on_parameters);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble();
    this->the_worker->move(next_ensemble,
                           current_ensemble,
                           conditioned_on_parameters);
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

MoveOutput* EnsembleKalmanFilter::move(RandomNumberGenerator &rng,
                                          Particle &ensemble_member,
                                          const Parameters &conditioned_on_parameters)
{
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    VectorSingleIndex index(this->current_index);
    if (i==0)
    {
      moved_ensemble_member = this->proposal_kernel->move(rng,
                                                          ensemble_member,
                                                          conditioned_on_parameters);
    }
    else
    {
      moved_ensemble_member = this->proposal_kernel->move(rng,
                                                          moved_ensemble_member,
                                                          conditioned_on_parameters);
    }
  }
  
  moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
  
  MoveOutput* moved_output = new SinglePointMoveOutput(moved_ensemble_member);
  
  return moved_output;
}

MoveOutput* EnsembleKalmanFilter::subsample_move(RandomNumberGenerator &rng,
                                                 Particle &ensemble_member,
                                                 const Parameters &conditioned_on_parameters)
{
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    VectorSingleIndex index(this->current_index);
    if (i==0)
    {
      moved_ensemble_member = this->proposal_kernel->subsample_move(rng,
                                                                    ensemble_member,
                                                                    conditioned_on_parameters);
    }
    else
    {
      moved_ensemble_member = this->proposal_kernel->subsample_move(rng,
                                                                    moved_ensemble_member,
                                                                    conditioned_on_parameters);
    }
  }
  
  moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
  
  MoveOutput* moved_output = new SinglePointMoveOutput(moved_ensemble_member);
  
  return moved_output;
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state,
                                                    const Parameters &conditioned_on_parameters)
{
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("EnsembleKalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->ensemble_kalman_evaluate_smcfixed_part(current_state,
                                               conditioned_on_parameters);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                 conditioned_on_parameters);
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                  const Parameters &conditioned_on_parameters)
{
  //this->the_worker->smcfixed_weight(current_state->all_ensembles.back());
  //current_state->initialise_next_step();
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
                                                                                    const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  //if (!this->sequencer_limit_is_fixed)
  //  this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    this->the_worker->shift(&current_state->back(),1.0);
    this->the_worker->unpack(&current_state->back());
    
    /*
     if (current_state->all_particles.size()==1)
     {
     
     this->the_worker->weight(current_state->back());
     }
     else
     {
     if (this->dynamic_proposal_is_evaluated)
     {
     this->the_worker->pf_weight(current_state->back(),
     *(current_state->all_particles.end()-2),
     this->proposal_kernel);
     }
     else
     {
     this->the_worker->pf_weight(current_state->back(),
     *(current_state->all_particles.end()-2),
     NULL);
     }
     }
     
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
     */
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_simulate(current_state,
                                   conditioned_on_parameters);
    
  }
}

void EnsembleKalmanFilter::subsample_ensemble_kalman_simulate(EnsembleKalmanOutput* current_state,
                                                              const Parameters &conditioned_on_parameters)
{
  if (current_state->all_ensembles.size()==0)
  {
    VectorSingleIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index,
                            conditioned_on_parameters);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble();
    this->the_worker->subsample_move(next_ensemble,
                                     current_ensemble,
                                     conditioned_on_parameters);
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void EnsembleKalmanFilter::subsample_ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state,
                                                              const Parameters &conditioned_on_parameters)
{
  if (!this->last_index_is_fixed)
  {
    throw std::runtime_error("EnsembleKalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->subsample_ensemble_kalman_evaluate_smcfixed_part(current_state,
                                                         conditioned_on_parameters);
  this->subsample_ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                           conditioned_on_parameters);
}

void EnsembleKalmanFilter::subsample_ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                            const Parameters &conditioned_on_parameters)
{
  //this->the_worker->smcfixed_weight(current_state->all_ensembles.back());
  //current_state->initialise_next_step();
}

void EnsembleKalmanFilter::subsample_ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
                                                                                              const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  //if (!this->sequencer_limit_is_fixed)
  //  this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    this->the_worker->shift(&current_state->back(),1.0);
    this->the_worker->unpack(&current_state->back());
    
    /*
     if (current_state->all_particles.size()==1)
     {
     
     this->the_worker->weight(current_state->back());
     }
     else
     {
     if (this->dynamic_proposal_is_evaluated)
     {
     this->the_worker->pf_weight(current_state->back(),
     *(current_state->all_particles.end()-2),
     this->proposal_kernel);
     }
     else
     {
     this->the_worker->pf_weight(current_state->back(),
     *(current_state->all_particles.end()-2),
     NULL);
     }
     }
     
     current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
     */
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->ensemble_factors->set_data(this->current_index);
    
    this->subsample_ensemble_kalman_simulate(current_state,
                                             conditioned_on_parameters);
    
  }
}

//void EnsembleKalmanFilter::evaluate(EnsembleKalmanOutput* current_state,
//                            const Parameters &conditioned_on_parameters)
//{
  // unsure
  // copied from KF??????
  
  /*
  
  // use conditioned_on_parameters to set the next index to stop on
  if (!this->last_index_is_fixed)
  {
    this->first_index = this->last_index+1;
    this->current_index = this->first_index;
    this->last_index = size_t(conditioned_on_parameters[this->index_name][0]);
    if (this->first_index>this->last_index)
      throw std::runtime_error("EnsembleKalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
  }
  
  // Set predictor and updater with parameters.
  this->predictor->conditioned_on_parameters = conditioned_on_parameters;
  this->updater->conditioned_on_parameters = conditioned_on_parameters;
  
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
  */
//}

//void EnsembleKalmanFilter::subsample_evaluate(EnsembleKalmanOutput* current_state,
//                                    const Parameters &conditioned_on_parameters)
//{
  // unsure
  // copied from KF??????
  
  /*
   
   // use conditioned_on_parameters to set the next index to stop on
   if (!this->last_index_is_fixed)
   {
   this->first_index = this->last_index+1;
   this->current_index = this->first_index;
   this->last_index = size_t(conditioned_on_parameters[this->index_name][0]);
   if (this->first_index>this->last_index)
   throw std::runtime_error("EnsembleKalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
   }
   
   // Set predictor and updater with parameters.
   this->predictor->conditioned_on_parameters = conditioned_on_parameters;
   this->updater->conditioned_on_parameters = conditioned_on_parameters;
   
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
   */
//}

bool EnsembleKalmanFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}


//double EnsembleKalmanFilter::evaluate(const Parameters &parameters)
//{
//  return this->func(parameters,*this->data);
//}

// void EnsembleKalmanFilter::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }

/*
void EnsembleKalmanFilter::weight_for_adapting_sequence(Ensemble &current_particles,
                                                        double incremental_temperature)
{
  // no adaptation so don't do anything
  //this->the_worker->smcadaptive_given_smcfixed_evaluate_target(current_particles);
}
*/
