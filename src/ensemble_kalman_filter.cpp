#include "ensemble_kalman_filter.h"
#include "caller_tempering.h"
#include "ensemble_kalman_output.h"
#include "factors.h"
#include "single_point_move_output.h"
#include "proposal_kernel.h"
#include "ensemble_factor_variables.h"
#include "ensemble_kalman_filter.h"
#include "ensemble_kalman_worker.h"
#include "ensemble_factors.h"
#include "vector_index.h"
#include "single_point_move_output.h"
#include "independent_proposal_kernel.h"
#include "measurement_covariance_estimator.h"
#include "hmm_ensemble_factors.h"
#include "sequential_ensemble_kalman_worker.h"
#include "positive_smc_criterion.h"

namespace ilike
{
EnsembleKalmanFilter::EnsembleKalmanFilter()
:EnsembleKalman()
{
  this->index.reset();
}

EnsembleKalmanFilter::EnsembleKalmanFilter(RandomNumberGenerator* rng_in,
                                           size_t* seed_in,
                                           Data* data_in,
                                           size_t lag_in,
                                           //const std::string &state_name,
                                           //const arma::colvec &prior_mean_in,
                                           //const arma::mat &prior_covariance_in,
                                           const std::string &index_name_in,
                                           const std::string &time_name_in,
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
  this->time_name = time_name_in;
  this->time_diff_name = time_diff_name_in;
  this->first_index = first_index_in;
  this->last_index = last_index_in;
  this->predictions_per_update = predictions_per_update_in;
  this->update_time_step = update_time_step_in;
  this->current_time = initial_time_in;
  this->initial_time = initial_time_in;
  this->current_index = this->first_index;
  
  this->proposal.reset(prior_in);
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
  
  this->index.reset(new VectorIndex(indices));
  
  this->transition_model_kernel = transition_model_in;
  
  this->ensemble_factors.reset(new HMMEnsembleFactors(transition_model_in,
                                                  measurement_covariance_estimators));
  
  if (parallel_in==TRUE)
  {
  }
  else
  {
    this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
  }

  // needed to store schedule parameters
  this->sequencer = EnsembleSequencer();
  
}

EnsembleKalmanFilter::~EnsembleKalmanFilter()
{
  
  this->index.reset();
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

  this->index.reset();
  
  EnsembleKalman::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* EnsembleKalmanFilter::duplicate() const
{
  return( new EnsembleKalmanFilter(*this));
}

EnsembleKalman* EnsembleKalmanFilter::ensemble_kalman_duplicate() const
{
  return( new EnsembleKalmanFilter(*this));
}

void EnsembleKalmanFilter::make_copy(const EnsembleKalmanFilter &another)
{
  //if (another.proposal_kernel!=NULL)
  // the transition kernel is owned by this filter's own HMMEnsembleFactors
  // (duplicated by EnsembleKalman::make_copy), not by the filter copied from
  this->transition_model_kernel = NULL;
  if (another.transition_model_kernel!=NULL)
    this->transition_model_kernel = dynamic_cast<HMMEnsembleFactors*>(this->ensemble_factors.get())->get_transition_kernel();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
  
  this->index_name = another.index_name;
  this->time_name = another.time_name;
  this->time_diff_name = another.time_diff_name;
  this->first_index = another.first_index;
  this->last_index = another.last_index;
  this->predictions_per_update = another.predictions_per_update;
  this->update_time_step = another.update_time_step;
  this->current_time = another.current_time;
  this->initial_time = another.initial_time;
  this->current_index = another.current_index;
  this->follows_caller_target = another.follows_caller_target;
  this->followed_variable = another.followed_variable;
  this->followed_anneal = another.followed_anneal;
  this->exchange_below = another.exchange_below;
  this->exchange_multiply = another.exchange_multiply;
  this->exchange_max = another.exchange_max;
  this->followed_steps = another.followed_steps; // (the same counter: shared)
  this->lag = another.lag;

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
  this->current_index = this->first_index;
  this->ensemble_factors->set_data(this->current_index);
  
  this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
  
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag,
                                                          this->transform,
                                                          this->results_name);
  
  return output;
}

void EnsembleKalmanFilter::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state)
{
  if (current_state->all_ensembles.size()==0)
  {
    VectorIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
    
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
  else
  {
    this->current_index = this->first_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  }
  
  this->ensemble_kalman_evaluate_smcfixed_part(current_state);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state);
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state)
{
}

MoveOutput* EnsembleKalmanFilter::move(RandomNumberGenerator &rng,
                                       Particle &ensemble_member)
{
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    
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

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->ensemble_factors->set_data(this->current_index);
    this->current_time = this->current_time + this->update_time_step;
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    
    this->ensemble_kalman_simulate(current_state);
    
  }
}

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_run(const Parameters &conditioned_on_parameters)
{
  EnsembleKalmanOutput* current_state = this->ensemble_kalman_initialise(conditioned_on_parameters);
  this->ensemble_kalman_simulate(current_state, conditioned_on_parameters);
  this->ensemble_kalman_evaluate(current_state, conditioned_on_parameters);
  
  return current_state;
}

EnsembleKalmanOutput* EnsembleKalmanFilter::specific_ensemble_kalman_initialise(const Parameters &parameters)
{
  this->current_index = this->first_index;
  this->ensemble_factors->set_data(this->current_index);
  
  this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
  
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this,
                                                          this->lag,
                                                          this->transform,
                                                          this->results_name);
  
  return output;
}

void EnsembleKalmanFilter::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state,
                                                    const Parameters &conditioned_on_parameters)
{
  // (C265) a following filter draws its initial ensemble when the output is
  // first evaluated, at the first time (follow_from_start)
  if (this->follows_caller_target)
    return;
  if (current_state->all_ensembles.size()==0)
  {
    // (C265) a new run starts at the first index and the initial time (both
    // were left where the previous run ended: the initial ensemble's index,
    // and a time that grew from run to run, as the particle filter's did
    // before C256)
    this->current_index = this->first_index;
    this->current_time = this->initial_time;
    VectorIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index,
                            conditioned_on_parameters);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
    
    this->the_worker->move(next_ensemble,
                           current_ensemble);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

// (C264) as move (one prediction per update), the moved member built in a
// dropped ensemble's member: the same random numbers (the transition's
// simulate_into draws as its simulate does) and the same member
MoveOutput* EnsembleKalmanFilter::move_reusing(RandomNumberGenerator &rng,
                                               Particle &ensemble_member,
                                               MoveOutput* spare)
{
  SinglePointMoveOutput* output = dynamic_cast<SinglePointMoveOutput*>(spare);
  if ((output == NULL) || (this->predictions_per_update != 1))
  {
    delete spare;
    return this->move(rng, ensemble_member);
  }
  this->transition_model_kernel->move_into(rng, ensemble_member, NULL, output->back());
  return output;
}

MoveOutput* EnsembleKalmanFilter::subsample_move(RandomNumberGenerator &rng,
                                                 Particle &ensemble_member)
{
  
  Particle moved_ensemble_member;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    
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

  MoveOutput* moved_output = new SinglePointMoveOutput(std::move(moved_ensemble_member));
  
  return moved_output;
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state,
                                                    const Parameters &conditioned_on_parameters)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("EnsembleKalmanFilter::evaluate - setting limit with parameters not yet written.");
  }
  else
  {
    this->current_index = this->first_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
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
  else
  {
    this->current_index = this->first_index;
    this->ensemble_factors->set_data(this->current_index);
    
    this->sequencer.schedule_parameters = Parameters();//this->time_name,this->current_time);
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
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
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();

    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->ensemble_factors->set_data(this->current_index);
    this->current_time = this->current_time + this->update_time_step;
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    
    this->ensemble_kalman_subsample_simulate(current_state,
                                             conditioned_on_parameters);
    
  }
  
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                  const Parameters &conditioned_on_parameters)
{
}

void EnsembleKalmanFilter::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
                                                                                    const Parameters &conditioned_on_parameters)
{
  // (C269, PMMH P1) the transition's values that depend only on its declared
  // parameters, for this run
  if (this->transition_model_kernel != NULL)
    this->transition_model_kernel->prepare_for_parameters(conditioned_on_parameters);
  if (this->follows_caller_target)
  {
    this->follow_evaluate(current_state, conditioned_on_parameters);
    return;
  }
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
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
    {
      current_state->write(results_name);
    }
    
    current_state->start_time = std::chrono::high_resolution_clock::now();

    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->current_index = this->current_index+1;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    
    this->ensemble_factors->set_data(this->current_index);
    this->current_time = this->current_time + this->update_time_step;
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    
    this->ensemble_kalman_simulate(current_state,
                                   conditioned_on_parameters);
    
  }
  
}

void EnsembleKalmanFilter::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* current_state)
{
  if (current_state->all_ensembles.size()==0)
  {
    VectorIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
    
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
    VectorIndex index(this->current_index);
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            &index,
                            conditioned_on_parameters);
  }
  else
  {
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
    
    this->the_worker->subsample_move(next_ensemble,
                                     current_ensemble);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

//void EnsembleKalmanFilter::evaluate(EnsembleKalmanOutput* current_state,
//                            const Parameters &conditioned_on_parameters)
//{
// unsure
// copied from KF??????

//void EnsembleKalmanFilter::subsample_evaluate(EnsembleKalmanOutput* current_state,
//                                    const Parameters &conditioned_on_parameters)
//{
// unsure
// copied from KF??????

bool EnsembleKalmanFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}

// --- SMC^2 with an ensemble Kalman filter (dev/plan.md 3.21 (c2''')) ----------

void EnsembleKalmanFilter::set_follows_caller_target(const std::string &fraction_variable, bool anneal)
{
  this->follows_caller_target = true;
  this->followed_variable = fraction_variable;
  this->followed_anneal = anneal;
  this->smcfixed_flag = false; // (the value changes along the caller's sequence)
}

// the filter's shared state for an output's step: the index, time and data of
// observation step + 1 (as the full run sets them)
void EnsembleKalmanFilter::follow_set_time_step(size_t step)
{
  this->current_index = this->first_index + step;
  this->current_time = this->initial_time + double(step)*this->update_time_step;
  this->ensemble_factors->set_data(this->current_index);
  this->sequencer.schedule_parameters = Parameters();
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
}

// the next observation's log-likelihood increment from the output's ensemble
// (already moved to that time; its measurement noise drawn with the move):
// packing and the measurement covariances are recomputed from the members,
// so a trial can be repeated
double EnsembleKalmanFilter::follow_increment(EnsembleKalmanOutput* current_state)
{
  this->the_worker->pack(&current_state->back());
  this->find_measurement_covariances(current_state);
  return current_state->calculate_latest_log_normalising_constant_ratio();
}

double EnsembleKalmanFilter::follow_trial_value(EnsembleKalmanOutput* current_state,
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
      this->follow_set_time_step(current_state->followed_steps);
      return value + ((r > 0.0) ? r : 1.0)*this->follow_increment(current_state);
    default:
      Rcpp::stop(caller_trial_too_far_message());
  }
  return value;
}

// keep the next observation: its increment, the update (shift) and unpack,
// then (unless it was the last) the move to the next time
void EnsembleKalmanFilter::follow_commit(EnsembleKalmanOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  size_t step = current_state->followed_steps;
  this->follow_set_time_step(step);
  current_state->committed_log_likelihood += this->follow_increment(current_state);
  this->the_worker->shift(&current_state->back(), 1.0);
  this->the_worker->unpack(&current_state->back());
  current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
  current_state->log_likelihood = current_state->committed_log_likelihood;
  current_state->set_llhd(current_state->committed_log_likelihood);
  current_state->set_time();
  current_state->followed_steps = step + 1;
  if (current_state->followed_steps > *this->followed_steps)
    *this->followed_steps = current_state->followed_steps;
  if ((current_state->followed_steps < this->observations()) && (current_state->committed_log_likelihood != -arma::datum::inf))
  {
    this->follow_set_time_step(step + 1);
    Ensemble* current_ensemble = &current_state->back();
    Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
    this->the_worker->move(next_ensemble, current_ensemble);
    current_state->increment_enk_iteration();
  }
}

// a new output: the initial ensemble at the first time, then the
// observations committed so far, then a trial at the caller's target
void EnsembleKalmanFilter::follow_from_start(EnsembleKalmanOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  size_t k;
  double r;
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  size_t committed = std::min(*this->followed_steps, k);
  current_state->followed_steps = 0;
  current_state->committed_log_likelihood = 0.0;
  current_state->log_likelihood = 0.0;
  this->follow_set_time_step(0);
  if (current_state->all_ensembles.size()==0)
  {
    VectorIndex index(this->current_index);
    this->simulate_proposal(current_state, &index, conditioned_on_parameters);
  }
  while ((current_state->followed_steps < committed) && (current_state->committed_log_likelihood != -arma::datum::inf))
    this->follow_commit(current_state, conditioned_on_parameters);
  current_state->caught_up = true;
  current_state->log_likelihood = this->follow_trial_value(current_state, conditioned_on_parameters, k, r);
}

void EnsembleKalmanFilter::follow_evaluate(EnsembleKalmanOutput* current_state,
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

// the caller's move stage: commit the next observation once it is whole; the
// value stays the one at the caller's current target
void EnsembleKalmanFilter::advance_at_caller_target_enk(EnsembleKalmanOutput* current_state,
                                                        const Parameters &conditioned_on_parameters)
{
  if (!this->follows_caller_target || !current_state->caught_up)
    return;
  if (this->transition_model_kernel != NULL)
    this->transition_model_kernel->prepare_for_parameters(conditioned_on_parameters); // (C269)
  size_t k;
  double r;
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  if ((current_state->committed_log_likelihood != -arma::datum::inf) && (k == current_state->followed_steps + 1) && (r == 0.0))
    this->follow_commit(current_state, conditioned_on_parameters);
  current_state->log_likelihood = this->follow_trial_value(current_state, conditioned_on_parameters, k, r);
}

// --- SMC^2's exchange step for the ensemble size (C267) ----------------------

void EnsembleKalmanFilter::set_exchange(double below, double multiply, size_t maximum)
{
  this->exchange_below = below;
  this->exchange_multiply = multiply;
  this->exchange_max = maximum;
}

size_t EnsembleKalmanFilter::exchange_size(double acceptance_rate, const Parameters &caller_parameters) const
{
  if (!this->has_exchange() || !this->follows_caller_target)
    return 0;
  size_t k;
  double r;
  caller_tempering_point(caller_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
  return caller_exchange_size(this->number_of_ensemble_members, acceptance_rate, this->exchange_below,
                              this->exchange_multiply, this->exchange_max, r);
}

// a new ensemble size for the outputs made from now on (each output copies
// this filter when it is made): the worker's per-member vector follows it
void EnsembleKalmanFilter::set_number_of_inner_particles(size_t number_of_inner_particles_in)
{
  if (number_of_inner_particles_in == this->number_of_ensemble_members)
    return;
  this->number_of_ensemble_members = number_of_inner_particles_in;
  SequentialEnsembleKalmanWorker* worker = dynamic_cast<SequentialEnsembleKalmanWorker*>(this->the_worker.get());
  if (worker == NULL)
    Rcpp::stop("SMC^2: the exchange step needs the ensemble filter's sequential worker.");
  worker->resize_for_members();
}

// (C270) this filter and its measurement estimators draw from rng_in, seed_in
void EnsembleKalmanFilter::redirect_random_numbers(RandomNumberGenerator* rng_in, size_t* seed_in)
{
  LikelihoodEstimator::redirect_random_numbers(rng_in, seed_in);
  if (this->ensemble_factors != NULL)
    this->ensemble_factors->redirect_random_numbers(rng_in, seed_in);
}

}
