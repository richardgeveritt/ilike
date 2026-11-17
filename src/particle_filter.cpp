#include <numeric>
#include "particle_filter.h"
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

void ParticleFilter::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
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

}
