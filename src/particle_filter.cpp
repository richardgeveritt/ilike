#include "vector_factors.h"
#include <numeric>
#include "particle_filter.h"
#include "caller_tempering.h"
#include "factors.h"
#include "chunk_parallel.h"
#include "progress_reporter.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "positive_smc_criterion.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "custom_distribution_proposal_kernel.h"
#include "hmm_factors.h"
#include "hmm_index.h"
#include "vector_index.h"

namespace ilike
{
ParticleFilter::ParticleFilter()
:SMC()
{
  this->index.reset();
}

ParticleFilter::ParticleFilter(RandomNumberGenerator* rng_in,
                               size_t* seed_in,
                               Data* data_in,
                               const Parameters &algorithm_parameters,
                               size_t number_of_particles_in,
                               size_t lag_in,
                               size_t lag_proposed_in,
                               const std::string &index_name_in,
                               const std::string &time_name_in,
                               const std::string &time_diff_name_in,
                               size_t first_index_in,
                               size_t last_index_in,
                               size_t predictions_per_update_in,
                               double update_time_step_in,
                               double initial_time_in,
                               SMCCriterion* adaptive_resampling_in,
                               const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                               IndependentProposalKernel* proposal_in,
                               ProposalKernel* transition_model_in,
                               ProposalKernel* transition_proposal_in,
                               VectorIndex* evaluated_in_initial_weight_update,
                               VectorIndex* evaluated_in_pf_weight_update,
                               bool proposal_is_evaluated_in,
                               bool transition_proposal_is_evaluated_in,
                               bool smcfixed_flag_in,
                               bool sequencer_limit_is_fixed_in,
                               bool transform_proposed_particles,
                               bool parallel_in,
                               size_t grain_size_in,
                               const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, algorithm_parameters, number_of_particles_in, lag_in, lag_proposed_in, transition_proposal_in->get_proposals(), adaptive_resampling_in, proposal_is_evaluated_in, smcfixed_flag_in, sequencer_limit_is_fixed_in, transform_proposed_particles, results_name_in)
{
  
  proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  transition_proposal_in->set_proposal_parameters(&this->algorithm_parameters);
  
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
  
  this->index.reset(new HMMIndex(evaluated_in_initial_weight_update->get_indices(),
                             evaluated_in_pf_weight_update->get_indices(),
                             transition_proposal_is_evaluated_in,
                             this->first_index));
  
  this->factors.reset(new HMMFactors(transition_model_in,
                                 likelihood_estimators_in));
  
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal_in,
                                                            likelihood_estimators_in));
  
  this->transition_proposal.reset(transition_proposal_in);
  
  this->parallel = parallel_in;
  this->chunk_size = grain_size_in;
  this->the_worker.reset(new SequentialSMCWorker(this));
  
  this->transition_proposal_is_evaluated = transition_proposal_is_evaluated_in;

  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new PositiveSMCCriterion();
  this->sequencer = Sequencer(this->the_worker.get(),
                              schedule_in,
                              "",
                              25,
                              smc_criterion);
  
  // also need to set up first proposal
  
}

//Copy constructor for the ParticleFilter class.
ParticleFilter::ParticleFilter(const ParticleFilter &another)
:SMC(another)
{
  this->make_copy(another);
}

//Destructor for the ParticleFilter class.
ParticleFilter::~ParticleFilter()
{
  this->transition_proposal.reset();
  
  this->index.reset();
}

void ParticleFilter::operator=(const ParticleFilter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->transition_proposal.reset();
  
  this->index.reset();
  
  SMC::operator=(another);
  this->make_copy(another);
}

SMC* ParticleFilter::smc_duplicate() const
{
  return( new ParticleFilter(*this));
}

LikelihoodEstimator* ParticleFilter::duplicate() const
{
  return( new ParticleFilter(*this));
}

void ParticleFilter::make_copy(const ParticleFilter &another)
{
  // B29: the copy owns its own transition proposal (the destructor deletes it)
  if (another.transition_proposal!=NULL)
    this->transition_proposal.reset(another.transition_proposal->proposal_kernel_duplicate());
  else
    this->transition_proposal.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->hmm_index_duplicate());
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
  this->follows_caller_target = another.follows_caller_target;
  this->followed_variable = another.followed_variable;
  this->followed_anneal = another.followed_anneal;
  this->exchange_below = another.exchange_below;
  this->exchange_multiply = another.exchange_multiply;
  this->exchange_max = another.exchange_max;
  this->followed_steps = another.followed_steps;
  this->current_index = another.current_index;
  this->lag = another.lag;
  this->transition_proposal_is_evaluated = another.transition_proposal_is_evaluated;
}

SMCOutput* ParticleFilter::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  ProgressScope progress(this->show_progress,
                         (this->last_index >= this->first_index) ? double(this->last_index - this->first_index + 1) : -1.0,
                         "Particle filter");
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* ParticleFilter::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->current_time = this->initial_time; // (it carried over from the previous run: c2')
  this->factors->set_data(this->current_index);
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->index->set_time_index(this->current_index);
  return output;
}

void ParticleFilter::simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update proposals
    this->transition_proposal->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state);
    //this->sequencer.find_next_target_bisection(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    this->increment_time_index();
    
    if (this->can_move_in_chunks())
    {
      this->move_in_chunks(next_particles,
                           current_particles);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->the_worker->move(next_particles,
                             current_particles);
    }
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::evaluate_smc(SMCOutput* current_state)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("ParticleFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  else
  {
    this->current_index = this->first_index;
    this->current_time = this->initial_time; // (it carried over from the previous run: c2')
    this->factors->set_data(this->current_index);
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->index->set_time_index(this->current_index);
  }
  
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void ParticleFilter::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
  
}

bool ParticleFilter::reuses_particle_objects() const
{
  return this->predictions_per_update == 1;
}

MoveOutput* ParticleFilter::move_reusing(RandomNumberGenerator &rng,
                                         const Particle &particle,
                                         MoveOutput* spare) const
{
  // as move (one prediction per update), building the new particle in a
  // dropped generation's object
  SinglePointMoveOutput* output = dynamic_cast<SinglePointMoveOutput*>(spare);
  if ((output == NULL) || (this->predictions_per_update != 1))
  {
    delete spare;
    return this->move(rng, particle);
  }
  this->transition_proposal->set_data(this->factors->get_current_data());
  return this->move_reusing_with_data_set(rng, particle, spare);
}

MoveOutput* ParticleFilter::move_reusing_with_data_set(RandomNumberGenerator &rng,
                                                       const Particle &particle,
                                                       MoveOutput* spare) const
{
  SinglePointMoveOutput* output = dynamic_cast<SinglePointMoveOutput*>(spare);
  if ((output == NULL) || (this->predictions_per_update != 1))
  {
    delete spare;
    return this->move_with_data_set(rng, particle);
  }
  this->transition_proposal->move_into(rng, particle, NULL, output->back());
  return output;
}

MoveOutput* ParticleFilter::move(RandomNumberGenerator &rng,
                                 const Particle &particle) const
{
  this->transition_proposal->set_data(this->factors->get_current_data());
  return this->move_with_data_set(rng, particle);
}

MoveOutput* ParticleFilter::move_with_data_set(RandomNumberGenerator &rng,
                                               const Particle &particle) const
{
  Particle moved_particle;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    if (i==0)
    {
      moved_particle = this->transition_proposal->move(rng,
                                                       particle);
    }
    else
    {
      moved_particle = this->transition_proposal->move(rng,
                                                       moved_particle);
    }
    
  }
  
  return new SinglePointMoveOutput(std::move(moved_particle));
}

void ParticleFilter::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    Rcpp::stop("ParticleFilter::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    if (this->show_progress)
      ProgressReporter::instance().update(double(this->current_index - this->first_index + 1));
    
    if (this->parallel && this->can_move_in_chunks())
    {
      this->weight_in_chunks(current_state);
    }
    else if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->weight(this->index.get(),
                               current_state->back());
    }
    else
    {
      
      if (this->transition_proposal_is_evaluated)
      {
        this->the_worker->pf_weight(this->index.get(),
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    this->transition_proposal.get());
      }
      else
      {
        this->the_worker->pf_weight(this->index.get(),
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // need to add Sequencer here when we add MCMC moves
    
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    current_state->set_time_and_reset_start();

    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state);
    
  }
}

//void ParticleFilter::weight_for_adapting_sequence(Particles &current_particles)
//{
// no adaptation so don't do anything
//this->the_worker->smcadaptive_given_smcfixed_evaluate_target(current_particles);
//}

SMCOutput* ParticleFilter::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* ParticleFilter::specific_initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->current_time = this->initial_time; // (it carried over from the previous run: c2')
  this->factors->set_data(this->current_index);
  this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->index->set_time_index(this->current_index);
  return output;
}

void ParticleFilter::simulate_smc(SMCOutput* current_state,
                                  const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update proposals
    this->transition_proposal->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    this->increment_time_index();
    
    if (this->can_move_in_chunks())
    {
      this->move_in_chunks(next_particles,
                           current_particles);
    }
    else
    {
      if (this->parallel)
        this->parallel_refused = true;
      this->the_worker->move(next_particles,
                             current_particles);
    }
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::evaluate_smc(SMCOutput* current_state,
                                  const Parameters &conditioned_on_parameters)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("ParticleFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  else
  {
    this->current_index = this->first_index;
    this->current_time = this->initial_time; // (it carried over from the previous run: c2')
    this->factors->set_data(this->current_index);
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->index->set_time_index(this->current_index);
  }
  
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void ParticleFilter::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
  
}

bool ParticleFilter::defers_simulation() const
{
  return !this->follows_caller_target && this->sequencer_limit_is_fixed &&
    !this->transition_proposal_is_evaluated && !this->proposal_is_evaluated;
}

// The likelihood estimate is the product over the steps of the mean weight,
// and in a bootstrap filter a step's weights are the measurement density at
// the particles: each mean is at most the density's maximum over the state.
// So after step t the final log estimate is at most the log estimate so far
// plus the logs of the maxima at the later steps.
double ParticleFilter::run_with_threshold(SMCOutput* simulation,
                                          const Parameters &conditioned_on_parameters,
                                          double threshold,
                                          bool &below)
{
  below = false;
  this->transition_proposal->prepare_for_parameters(conditioned_on_parameters); // (C269)
  this->simulate_smc(simulation, conditioned_on_parameters);
  if (!this->defers_simulation() || (this->factors == NULL) || (this->last_index < this->first_index))
    return simulation->evaluate(conditioned_on_parameters);
  // the maxima at each step (the factors' data and index set for that step)
  size_t T = this->last_index - this->first_index + 1;
  this->stop_remaining.assign(T + 1, 0.0);
  Parameters at_step = conditioned_on_parameters.deep_copy();
  for (size_t k = 0; k < T; ++k)
  {
    size_t step = this->first_index + k;
    this->factors->set_data(step);
    this->index->set_time_index(step);
    if (!this->index_name.empty())
      at_step[this->index_name] = arma::mat(1, 1, arma::fill::value(double(step)));
    this->stop_remaining[k] = this->factors->log_maximum(this->index.get(), at_step);
  }
  for (size_t k = T; k > 0; --k)
    this->stop_remaining[k - 1] = this->stop_remaining[k - 1] + this->stop_remaining[k];
  this->stop_threshold = threshold;
  this->threshold_stopped = false;
  this->threshold_active = true;
  double value = simulation->evaluate(conditioned_on_parameters); // (sets the data and index back to the first step)
  this->threshold_active = false;
  if (this->threshold_stopped)
  {
    below = true;
    return -arma::datum::inf;
  }
  return value;
}

void ParticleFilter::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                  const Parameters &conditioned_on_parameters)
{
  // (C269, PMMH P1) the transition's values that depend only on its declared
  // parameters, for this run (a sequential point: before any parallel chunk)
  this->transition_proposal->prepare_for_parameters(conditioned_on_parameters);
  if (this->follows_caller_target)
  {
    // (c2') a new output is run through the observations committed so far;
    // afterwards an evaluation is a trial at the caller's current target:
    // the committed observations, plus the next one whole or (F1) a share
    if (!current_state->caught_up)
    {
      this->follow_from_start(current_state, conditioned_on_parameters);
      return;
    }
    size_t k;
    double r;
    this->tempering_point(conditioned_on_parameters, k, r);
    double value = this->trial_value(current_state, k, r);
    current_state->log_likelihood = value;
    current_state->set_latest_llhd(value);
    return;
  }
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    if (this->show_progress)
      ProgressReporter::instance().update(double(this->current_index - this->first_index + 1));
    
    if (this->parallel && this->can_move_in_chunks())
    {
      this->weight_in_chunks(current_state);
    }
    else if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->weight(this->index.get(),
                               current_state->back());
    }
    else
    {
      
      if (this->transition_proposal_is_evaluated)
      {
        this->the_worker->pf_weight(this->index.get(),
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    this->transition_proposal.get());
        
      }
      else
      {
        this->the_worker->pf_weight(this->index.get(),
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

    current_state->set_time_and_reset_start();

    // (dev/plan.md 3.21) stop when the estimate cannot reach the threshold:
    // the log estimate so far plus the maxima of the steps not yet weighted
    if (this->threshold_active)
    {
      size_t done = this->current_index - this->first_index + 1;
      double remaining = (done < this->stop_remaining.size()) ? this->stop_remaining[done] : 0.0;
      if (current_state->log_likelihood + remaining < this->stop_threshold)
      {
        this->threshold_stopped = true;
        current_state->terminate();
        break;
      }
    }

    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

void ParticleFilter::subsample_simulate_smc(SMCOutput* current_state,
                                            const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update proposals
    this->transition_proposal->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.subsample_find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    this->increment_time_index();
    
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::subsample_evaluate_smc(SMCOutput* current_state)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("ParticleFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  else
  {
    this->current_index = this->first_index;
    this->current_time = this->initial_time; // (it carried over from the previous run: c2')
    this->factors->set_data(this->current_index);
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->index->set_time_index(this->current_index);
  }
  
  this->subsample_evaluate_smcfixed_part_smc(current_state);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void ParticleFilter::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void ParticleFilter::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    if (this->show_progress)
      ProgressReporter::instance().update(double(this->current_index - this->first_index + 1));
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->subsample_weight(this->index.get(),
                                         current_state->back());
      
    }
    else
    {
      if (this->transition_proposal_is_evaluated)
      {
        this->the_worker->subsample_pf_weight(this->index.get(),
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              this->transition_proposal.get());
      }
      else
      {
        this->the_worker->subsample_pf_weight(this->index.get(),
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state);
    
  }
}

void ParticleFilter::subsample_simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update proposals
    this->transition_proposal->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.subsample_find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    this->increment_time_index();
    
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    this->subsample_evaluate_smcfixed_part_smc(current_state);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::subsample_evaluate_smc(SMCOutput* current_state,
                                            const Parameters &conditioned_on_parameters)
{
  if (!this->sequencer_limit_is_fixed)
  {
    Rcpp::stop("ParticleFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  else
  {
    this->current_index = this->first_index;
    this->current_time = this->initial_time; // (it carried over from the previous run: c2')
    this->factors->set_data(this->current_index);
    this->sequencer.schedule_parameters[this->time_diff_name] = this->update_time_step/double(this->predictions_per_update);
    this->sequencer.schedule_parameters[this->time_name] = this->current_time;
    this->sequencer.schedule_parameters[this->index_name] = this->current_index;
    this->index->set_time_index(this->current_index);
  }
  
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void ParticleFilter::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                          const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void ParticleFilter::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                            const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  else
    this->sequencer.reset();
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    check_interrupt(); // between SMC steps (dev/plan.md 3.11)
    if (this->show_progress)
      ProgressReporter::instance().update(double(this->current_index - this->first_index + 1));
    
    if (current_state->number_of_smc_iterations()==1)
    {
      
      this->the_worker->subsample_weight(this->index.get(),
                                         current_state->back());

    }
    else
    {
      if (this->transition_proposal_is_evaluated)
      {
        this->the_worker->subsample_pf_weight(this->index.get(),
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              this->transition_proposal.get());
      }
      else
      {
        this->the_worker->subsample_pf_weight(this->index.get(),
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

void ParticleFilter::weight_for_adapting_sequence(const Index* index,
                                                  Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_evaluate_target(index,
                                                               current_particles);
}

MoveOutput* ParticleFilter::subsample_move(RandomNumberGenerator &rng,
                                           const Particle &particle) const
{
  
  Particle moved_particle;
  this->transition_proposal->set_data(this->factors->get_current_data());
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    if (i==0)
    {
      moved_particle = this->transition_proposal->subsample_move(rng,
                                                                 particle);
    }
    else
    {
      moved_particle = this->transition_proposal->subsample_move(rng,
                                                                 moved_particle);
    }
  }
  return new SinglePointMoveOutput(std::move(moved_particle));
}

void ParticleFilter::subsample_weight_for_adapting_sequence(const Index* index,
                                                            Particles &current_particles)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_evaluate_target(index,
                                                                         current_particles);
}

bool ParticleFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}

void ParticleFilter::increment_time_index()
{
  // should move to sequencer
  this->current_index = this->current_index+1;
  
  // should move to sequencer
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  
  this->factors->set_data(this->current_index);
  this->index->set_time_index(this->current_index);
  this->current_time = this->current_time + this->update_time_step;
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
}

bool ParticleFilter::supports_parallel_copies() const
{
  return (this->compact_stored_output || (this->lag == 0)) &&
    (this->factors!=NULL) && this->factors->supports_parallel_copies();
}

bool ParticleFilter::can_move_in_chunks() const
{
  return (this->factors!=NULL) && this->factors->supports_parallel_chunks() &&
    (this->transition_proposal!=NULL) &&
    (dynamic_cast<const SequentialSMCWorker*>(this->the_worker.get())!=NULL);
}

void ParticleFilter::move_in_chunks(Particles* next_particles,
                                    Particles* current_particles)
{
  // SequentialSMCWorker::specific_move and the seed bookkeeping of
  // SMCWorker::move, grain_size particles at a time: each chunk with its own
  // stream, writing only its own particles' slots (spares, the dropped
  // generation's objects, are handed out by index).
  SMCWorker* worker = this->the_worker.get();
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  this->transition_proposal->set_data(this->factors->get_current_data());

  std::vector<std::unique_ptr<MoveOutput>> spares;
  spares.swap(next_particles->spare);
  if (spares.size() < N)
    spares.resize(N);
  next_particles->particles.clear();
  next_particles->particles.resize(N);

  RandomNumberGenerator first_stream(*worker->get_rng());
  first_stream.seed(worker->get_seed());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first_stream, K);

  run_chunk_groups(0, K, this->parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     for (size_t k = first_chunk; k < end_chunk; ++k)
                     {
                       size_t first = k*B;
                       size_t n = std::min(B, N - first);
                       RandomNumberGenerator local_rng(streams[k]);
                       for (size_t i = first; i < first + n; ++i)
                       {
                         const Particle &ancestor = (*current_particles)[current_particles->ancestor_variables[i]]->back();
                         MoveOutput* spare = spares[i].release();
                         MoveOutput* moved = (spare != NULL) ?
                           this->move_reusing_with_data_set(local_rng, ancestor, spare) :
                           this->move_with_data_set(local_rng, ancestor);
                         moved->back().previous_self = &ancestor;
                         next_particles->particles[i].reset(moved);
                       }
                     }
                   });

  worker->set_seed(worker->get_seed() + N);
  this->rng->seed(worker->get_seed());
  next_particles->simulate_resampling_variables(*this->rng);
  worker->set_seed(worker->get_seed() + N);
}

void ParticleFilter::weight_in_chunks(SMCOutput* current_state)
{
  // SequentialSMCWorker::weight (first step) or pf_weight, by chunks of
  // particles; no random numbers, so the same values in any order
  size_t N = this->number_of_particles;
  size_t B = ((this->chunk_size == 0) || (this->chunk_size > N)) ? N : this->chunk_size;
  size_t K = (N + B - 1)/B;
  Particles &current = current_state->back();
  bool first_step = (current_state->number_of_smc_iterations()==1);
  Particles* previous = first_step ? NULL : &*(current_state->end()-2);
  const ProposalKernel* kernel = this->transition_proposal_is_evaluated ? this->transition_proposal.get() : NULL;
  arma::colvec weights(N);
  run_chunk_groups(0, K, this->parallel, this->number_of_threads,
                   [&](size_t group, size_t first_chunk, size_t end_chunk)
                   {
                     for (size_t i = first_chunk*B; i < std::min(N, end_chunk*B); ++i)
                     {
                       Particle &particle = current[i]->back();
                       if (first_step)
                       {
                         weights[i] = particle.evaluate_likelihoods(this->index.get()) - particle.previous_target_evaluated;
                       }
                       else
                       {
                         weights[i] = particle.evaluate_likelihoods(this->index.get());
                         if (kernel != NULL)
                           weights[i] = weights[i] - kernel->evaluate_kernel(particle,
                                                                            (*previous)[previous->ancestor_variables[i]]->back());
                       }
                     }
                   });
  this->the_worker->set_unnormalised_log_incremental_weights(weights);
}


// --- SMC^2: the filter following data tempering (dev/plan.md 3.21 (c2')) ----

void ParticleFilter::set_follows_caller_target(const std::string &fraction_variable, bool anneal)
{
  this->follows_caller_target = true;
  this->followed_variable = fraction_variable;
  this->followed_anneal = anneal;
  this->smcfixed_flag = false;
  this->force_rich_generations = true;
  this->compact_stored_output = false;
}

// --- SMC^2's exchange step (dev/plan.md 3.21 (c2''), F2) ---------------------

void ParticleFilter::set_exchange(double below, double multiply, size_t maximum)
{
  this->exchange_below = below;
  this->exchange_multiply = multiply;
  this->exchange_max = maximum;
}

// the new number of particles, or 0 for no exchange: only when the caller's
// acceptance rate is below the threshold, the number is below the limit, and
// the caller's target is at a whole observation (between two, its targets
// for different numbers of particles have different marginals: C258)
size_t ParticleFilter::exchange_size(double acceptance_rate, const Parameters &caller_parameters) const
{
  if (!this->has_exchange() || !this->follows_caller_target)
    return 0;
  size_t k;
  double r;
  this->tempering_point(caller_parameters, k, r);
  return caller_exchange_size(this->number_of_particles, acceptance_rate, this->exchange_below,
                              this->exchange_multiply, this->exchange_max, r);
}

// a new number of particles for the outputs made from now on: the worker's
// weights, the default ancestors, and the resampling threshold (a number of
// particles: ess(p) is p N) follow it
void ParticleFilter::set_number_of_particles(size_t number_of_particles_in)
{
  if (number_of_particles_in == this->number_of_particles)
    return;
  double ratio = double(number_of_particles_in)/double(this->number_of_particles);
  this->number_of_particles = number_of_particles_in;
  this->setup_default_ancestor_variables();
  SequentialSMCWorker* worker = dynamic_cast<SequentialSMCWorker*>(this->the_worker.get());
  if (worker == NULL)
    Rcpp::stop("SMC^2: the exchange step needs the filter's sequential worker.");
  worker->resize_for_particles();
  if (this->resampling_criterion)
    this->resampling_criterion->scale_desired_criterion(ratio);
}

// the caller's target: n fraction = k + r, k whole observations and a share
// r in [0, 1) of the next. Without anneal, k = ceil(n fraction) and r = 0, as
// data(…, anneal = false) adds whole data points to factors of model
// statements; with anneal, k = floor(n fraction), as plates do
// (tempering_point). A fraction written as j/n gives k = j, r = 0 in both
// (rounding allowed for)
void ParticleFilter::tempering_point(const Parameters &conditioned_on_parameters, size_t &k, double &r) const
{
  caller_tempering_point(conditioned_on_parameters[this->followed_variable][0], this->observations(),
                         this->followed_anneal, k, r);
}

// the value at the caller's target (k, r): the committed observations, plus
// the next one's weighting (the share r of it with anneal); no random numbers
double ParticleFilter::trial_value(SMCOutput *current_state, size_t k, double r)
{
  double value = current_state->committed_log_likelihood;
  if (value == -arma::datum::inf)
    return value;
  size_t committed = current_state->followed_steps;
  switch (caller_trial(committed, k, r))
  {
    case CallerTrial::committed:
      return value;
    case CallerTrial::next:
      this->set_time_step(this->first_index + committed);
      return value + this->weight_current(current_state, (r > 0.0) ? r : 1.0);
    default:
      Rcpp::stop(caller_trial_too_far_message());
  }
  return value;
}

void ParticleFilter::set_time_step(size_t step)
{
  this->current_index = step;
  this->current_time = this->initial_time + double(step - this->first_index)*this->update_time_step;
  this->sequencer.schedule_parameters[this->index_name] = this->current_index;
  this->sequencer.schedule_parameters[this->time_name] = this->current_time;
  this->factors->set_data(this->current_index);
  this->index->set_time_index(this->current_index);
}

// the current particles weighted by the observation at the current time (the
// weights start again from the step's starting weights: no state is kept);
// the log of the estimated ratio
//
// (F1) power r < 1: the log ratio times r. The caller's intermediate target
// is then the geometric bridge between the two extended targets, the start
// (the committed observations, with the particles already propagated: the
// same psi) and the end (one more observation): p(theta) psi Zhat_start^(1-r)
// Zhat_end^r, as between any two targets of a main sequence (geometric
// annealing in SMC with transformations, arXiv:1612.06468 §2.4.1). It needs
// nothing of the filter beyond its two estimates, so it holds whatever the
// filter's proposals are
double ParticleFilter::weight_current(SMCOutput *current_state, double power)
{
  if (current_state->number_of_smc_iterations() == 1)
    this->the_worker->weight(this->index.get(), current_state->back());
  else
    this->the_worker->pf_weight(this->index.get(), current_state->back(), *(current_state->end() - 2),
                                this->transition_proposal_is_evaluated ? this->transition_proposal.get() : NULL);
  if (power != 1.0)
  {
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    return power*current_state->calculate_latest_log_normalising_constant_ratio();
  }
  current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  return current_state->calculate_latest_log_normalising_constant_ratio();
}

// the generation's entries (one log likelihood and one time each)
static void record_filter_generation(SMCOutput *current_state)
{
  current_state->set_llhd(current_state->committed_log_likelihood);
  current_state->set_time_and_reset_start();
}

// keep the weighting by the next observation, then (unless it was the last)
// resample and propagate to the next time
void ParticleFilter::commit_observation(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  size_t k = current_state->followed_steps; // (committed so far)
  this->set_time_step(this->first_index + k);
  current_state->committed_log_likelihood += this->weight_current(current_state);
  current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
  current_state->set_latest_llhd(current_state->committed_log_likelihood);
  current_state->followed_steps = k + 1;
  if (current_state->followed_steps > this->followed_steps)
    this->followed_steps = current_state->followed_steps;
  if ((current_state->followed_steps < this->observations()) && (current_state->committed_log_likelihood != -arma::datum::inf))
  {
    this->simulate_smc(current_state, conditioned_on_parameters); // (resample, then increment_time_index and move)
    record_filter_generation(current_state);
  }
}

void ParticleFilter::advance_at_caller_target_smc(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  if (!this->follows_caller_target || !current_state->caught_up)
    return;
  this->transition_proposal->prepare_for_parameters(conditioned_on_parameters); // (C269)
  // (F1) an observation is committed only once whole: at a share r > 0 the
  // filter is left as it is (the trial weighting is recomputed at the next
  // reweighting)
  size_t k;
  double r;
  this->tempering_point(conditioned_on_parameters, k, r);
  if ((current_state->committed_log_likelihood != -arma::datum::inf) && (k == current_state->followed_steps + 1) && (r == 0.0))
    this->commit_observation(current_state, conditioned_on_parameters);
  if (r > 0.0)
  {
    // the value stays the one at the caller's current target, committed +
    // r log(ratio): the next reweighting's increments are measured from it
    // (recomputed, since the caller's last trial may have been at another
    // fraction; no random numbers)
    double value = this->trial_value(current_state, k, r);
    current_state->set_latest_llhd(value);
    current_state->log_likelihood = value;
  }
  else
    current_state->log_likelihood = current_state->committed_log_likelihood;
}

// a new output (its particles drawn at the first time): committed and moved
// through the observations committed so far (the same steps as an output
// that has followed from the caller's start: the PMMH ratio needs the same
// law), then a trial at the caller's current observation if it is one more
void ParticleFilter::follow_from_start(SMCOutput *current_state, const Parameters &conditioned_on_parameters)
{
  size_t k;
  double r;
  this->tempering_point(conditioned_on_parameters, k, r);
  size_t committed = std::min(this->followed_steps, k);
  current_state->followed_steps = 0;
  current_state->committed_log_likelihood = 0.0;
  record_filter_generation(current_state);
  while ((current_state->followed_steps < committed) && (current_state->committed_log_likelihood != -arma::datum::inf))
    this->commit_observation(current_state, conditioned_on_parameters);
  current_state->caught_up = true;
  bool trial = (current_state->committed_log_likelihood != -arma::datum::inf) &&
               ((k > current_state->followed_steps) || (r > 0.0));
  double value = this->trial_value(current_state, k, r);
  if (trial)
    current_state->set_latest_llhd(value);
  current_state->log_likelihood = value;
}

}
