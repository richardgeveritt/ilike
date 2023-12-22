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
#include "independent_proposal_kernel.h"
#include "measurement_covariance_estimator.h"
#include "hmm_ensemble_factors.h"
#include "sequential_ensemble_kalman_worker.h"
#include "positive_smc_criterion.h"
//#include "ensemble_kalman_updater.h"
//#include "ensemble_kalman_predictor.h"

EnsembleKalmanFilter::EnsembleKalmanFilter()
  :EnsembleKalman()
{
  //this->proposal_kernel = NULL;
  this->index = NULL;
}

EnsembleKalmanFilter::EnsembleKalmanFilter(RandomNumberGenerator* rng_in,
                                           size_t* seed_in,
                                           Data* data_in,
                                           size_t lag_in,
                                           //const std::string &state_name,
                                           //const arma::colvec &prior_mean_in,
                                           //const arma::mat &prior_covariance_in,
                                           const std::string &index_name_in,
                                           //const std::string &time_name_in,
                                           const std::string &time_diff_name_in,
                                           //const std::vector<std::string> &measurements_names_in,
                                           size_t first_index_in,
                                           size_t last_index_in,
                                           size_t predictions_per_update_in,
                                           double update_time_step_in,
                                           double initial_time_in,
                                           size_t number_of_ensemble_members_in,
                                           EnsembleShifter* shifter_in,
                                           std::shared_ptr<Transform> transform_in,
                                           IndependentProposalKernel* prior_in,
                                           ProposalKernel* transition_model_in,
                                           const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in,
                                           bool smcfixed_flag_in,
                                           bool sequencer_limit_is_fixed_in,
                                           bool parallel_in,
                                           size_t grain_size_in,
                                           const std::string &results_name_in)
:EnsembleKalman(rng_in,
                seed_in,
                data_in,
                number_of_ensemble_members_in,
                lag_in,
                shifter_in,
                transform_in,
                smcfixed_flag_in,
                sequencer_limit_is_fixed_in,
                results_name_in)
{
  // Go through ENKI stuff.

  this->index_name = index_name_in;
  //this->time_name = time_name_in;
  this->time_diff_name = time_diff_name_in;
  this->first_index = first_index_in;
  this->last_index = last_index_in;
  this->predictions_per_update = predictions_per_update_in;
  this->update_time_step = update_time_step_in;
  this->current_time = initial_time_in;
  this->current_index = this->first_index;
  
  this->proposal = prior_in;
  this->proposal->set_proposal_parameters(&this->algorithm_parameters);
  
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators = measurement_covariance_estimators_in;
  
  std::vector<size_t> indices;
  
  for (size_t i=0;
       i<measurement_covariance_estimators.size();
       ++i)
  {
    measurement_covariance_estimators[i]->change_data();
    indices.push_back(i);
  }
  
  this->index = new VectorSingleIndex(indices);
  
  this->transition_model_kernel = transition_model_in;
  
  this->ensemble_factors = new HMMEnsembleFactors(transition_model_in,
                                                  measurement_covariance_estimators);
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,
    //this->model_and_algorithm.particle_simulator,
    //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialEnsembleKalmanWorker(this);
  }
  
  /*
  std::vector<double> schedule_in;
  schedule_in.reserve(this->last_index-this->first_index+1);
  for (size_t i=0; i<this->last_index-this->first_index+1; ++i)
  {
    schedule_in.push_back(this->first_index+i);
  }
  
  this->sequencer = EnsembleSequencer(this->the_worker,
                                      schedule_in,
                                      this->index_name,
                                      100,
                                      new PositiveSMCCriterion());
  */
  
  // needed to store schedule parameters
  this->sequencer = EnsembleSequencer();
  
}

EnsembleKalmanFilter::~EnsembleKalmanFilter()
{
  //if (this->proposal_kernel!=NULL)
  //  delete this->proposal_kernel;
  
  if (this->index!=NULL)
    delete this->index;
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
  
  if (this->index!=NULL)
    delete this->index;
  
  EnsembleKalman::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* EnsembleKalmanFilter::duplicate() const
{
  return( new EnsembleKalmanFilter(*this));
}

void EnsembleKalmanFilter::make_copy(const EnsembleKalmanFilter &another)
{
  //if (another.proposal_kernel!=NULL)
  this->transition_model_kernel = another.transition_model_kernel;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
  
  this->index_name = another.index_name;
  //this->measurements_names = another.measurements_names;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->current_time = another.current_time;
  this->current_index = another.current_index;
  //this->last_index_is_fixed = another.last_index_is_fixed;
  this->lag = another.lag;
  
  //this->schedule_parameters = another.schedule_parameters;
  
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

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_ensemble_kalman_initialise()
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag,
                                                          this->transform,
                                                          this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->ensemble_factors->set_data(this->current_index);
  
  this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  
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
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
    
    this->the_worker->move(next_ensemble,
                           current_ensemble);
    
    current_state->increment_enk_iteration();
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("EnsembleKalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
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
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    this->current_time = this->current_time + predict_time_step;
    
    VectorSingleIndex index(this->current_index);
    if (i==0)
    {
      moved_ensemble_member = this->transition_model_kernel->move(rng,
                                                                     ensemble_member);
    }
    else
    {
      moved_ensemble_member = this->transition_model_kernel->move(rng,
                                                                     moved_ensemble_member);
    }
  }
  
  //moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
    
  MoveOutput* moved_output = new SinglePointMoveOutput(std::move(moved_ensemble_member));
  
  return moved_output;
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    Rcpp::stop("EnsembleKalmanFilter::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    // will need to change sequencer to have temperature within t - at the moment temperature will not be set correctly
    
    this->the_worker->shift(&current_state->back(),1.0);
    this->the_worker->unpack(&current_state->back());
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->llhds.push_back(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_simulate(current_state);
    
  }
}

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_run(const Parameters &conditioned_on_parameters)
{
  EnsembleKalmanOutput* current_state = this->ensemble_kalman_initialise(conditioned_on_parameters);
  this->ensemble_kalman_simulate(current_state, conditioned_on_parameters);
  this->ensemble_kalman_evaluate(current_state, conditioned_on_parameters);
  
  //EnsembleKalmanOutput* current_state = this->ensemble_kalman_initialise(conditioned_on_parameters);
  //this->ensemble_kalman_evaluate(current_state, conditioned_on_parameters);
  //this->ensemble_kalman_evaluate(current_state);
  return current_state;
}

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_ensemble_kalman_initialise(const Parameters &parameters)
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag,
                                                          this->transform,
                                                          this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->factors->set_data(this->current_index);
  
  this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  
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
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
    
    /*
    this->the_worker->move(next_ensemble,
                           current_ensemble,
                           conditioned_on_parameters);
    */
    
    this->the_worker->move(next_ensemble,
                           current_ensemble);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

/*
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
  
  //moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
  
  MoveOutput* moved_output = new SinglePointMoveOutput(std::move(moved_ensemble_member));
  
  return moved_output;
}
*/

MoveOutput* EnsembleKalmanFilter::subsample_move(RandomNumberGenerator &rng,
                                                 Particle &ensemble_member)
{
  double predict_time_step = this->update_time_step/double(this->predictions_per_update);
  
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    this->current_time = this->current_time + predict_time_step;
    
    VectorSingleIndex index(this->current_index);
    if (i==0)
    {
      moved_ensemble_member = this->transition_model_kernel->subsample_move(rng,
                                                                               ensemble_member);
    }
    else
    {
      moved_ensemble_member = this->transition_model_kernel->subsample_move(rng,
                                                                               moved_ensemble_member);
    }
  }
  
  //moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
  
  MoveOutput* moved_output = new SinglePointMoveOutput(std::move(moved_ensemble_member));
  
  return moved_output;
}

/*
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
  
  //moved_ensemble_member.simulate_ensemble_factor_variables(&ensemble_member);
  
  MoveOutput* moved_output = new SinglePointMoveOutput(std::move(moved_ensemble_member));
  
  return moved_output;
}
*/

void EnsembleKalmanFilter::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state,
                                                    const Parameters &conditioned_on_parameters)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("EnsembleKalmanFilter::evaluate - setting limit with parameters not yet written.");
  }
  
  this->ensemble_kalman_evaluate_smcfixed_part(current_state,
                                               conditioned_on_parameters);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                 conditioned_on_parameters);
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* current_state,
                                                              const Parameters &conditioned_on_parameters)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("EnsembleKalmanFilter::evaluate - setting limit with parameters not yet written.");
  }
  
  this->ensemble_kalman_subsample_evaluate_smcfixed_part(current_state,
                                                         conditioned_on_parameters);
  this->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                           conditioned_on_parameters);
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                                            const Parameters &conditioned_on_parameters)
{
  
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
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
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    this->the_worker->shift(&current_state->back(),
                            1.0);
    this->the_worker->unpack(&current_state->back());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->llhds.push_back(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_subsample_simulate(current_state,
                                             conditioned_on_parameters);
    
  }

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
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    this->the_worker->shift(&current_state->back(),
                            1.0);
    this->the_worker->unpack(&current_state->back());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->llhds.push_back(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_simulate(current_state,
                                   conditioned_on_parameters);
    
  }

}

void EnsembleKalmanFilter::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* current_state)
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
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
    this->the_worker->subsample_move(next_ensemble,
                                     current_ensemble);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* current_state,
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
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
    //this->the_worker->subsample_move(next_ensemble,
    //                                 current_ensemble,
    //                                 conditioned_on_parameters);
    this->the_worker->subsample_move(next_ensemble,
                                     current_ensemble);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

/*
void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* current_state,
                                                              const Parameters &conditioned_on_parameters)
{
  if (!this->last_index_is_fixed)
  {
    Rcpp::stop("EnsembleKalmanFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->ensemble_kalman_subsample_evaluate_smcfixed_part(current_state,
                                                         conditioned_on_parameters);
  this->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                           conditioned_on_parameters);
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                            const Parameters &conditioned_on_parameters)
{
  //this->the_worker->smcfixed_weight(current_state->all_ensembles.back());
  //current_state->initialise_next_step();
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
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
    this->the_worker->shift(&current_state->back());
    this->the_worker->unpack(&current_state->back());
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->ensemble_factors->set_data(this->current_index);
    
    this->ensemble_kalman_subsample_simulate(current_state,
                                             conditioned_on_parameters);
    
  }
}
*/

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
      Rcpp::stop("EnsembleKalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
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
   Rcpp::stop("EnsembleKalmanFilter::evaluate - last index from parameters is before the current state of the filter.");
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

/*
void EnsembleKalmanFilter::setup_variables()
{
  this->setup_variables_using_candidate_parameters(this->proposal->independent_simulate(*this->rng));
}
*/

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
