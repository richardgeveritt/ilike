#include <numeric>
#include "particle_filter.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "positive_smc_criterion.h"
#include "move_output.h"
#include "single_point_move_output.h"
//#include "custom_independent_proposal_kernel.h"
#include "custom_distribution_proposal_kernel.h"
#include "hmm_factors.h"
#include "vector_single_index.h"

ParticleFilter::ParticleFilter()
   :SMC()
{

}

ParticleFilter::ParticleFilter(RandomNumberGenerator* rng_in,
                               size_t* seed_in,
                               Data* data_in,
                               size_t number_of_particles_in,
                               size_t lag_in,
                               size_t lag_proposed_in,
                               double resampling_desired_ess_in,
                               ProposalKernel* proposal_kernel_in,
                               EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                               EvaluateLogDistributionPtr evaluate_log_prior_in,
                               SimulateDistributionPtr simulate_proposal_in,
                               EvaluateLogDistributionPtr evaluate_log_proposal_in,
                               bool parallel_in,
                               size_t grain_size_in,
                               const std::string &results_name_in)
  :SMC(rng_in, seed_in, data_in, Parameters(), number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, resampling_desired_ess_in, true, false, true, results_name_in)
{
  
  proposal_kernel_in->set_proposal_parameters(&this->algorithm_parameters);
  
  // also need to set up first proposal
  
  /*
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  TRUE);
  
  EvaluateLogDistributionPtr power = annealing_power;
  
  this->model_and_algorithm.likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                                            seed_in,
                                                                                            data_in,
                                                                                            prior_times_likelihood,
                                                                                            power,
                                                                                            TRUE));
  
  // Proposal.
  ExactLikelihoodEstimator* proposal = new ExactLikelihoodEstimator(rng_in,
                                                                    seed_in,
                                                                    data_in,
                                                                    evaluate_log_proposal_in,
                                                                    TRUE);
  
  EvaluateLogDistributionPtr second_power = annealing_one_minus_power;
  
  this->model_and_algorithm.likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                                            seed_in,
                                                                                            data_in,
                                                                                            proposal,
                                                                                            second_power,
                                                                                            TRUE));
  
  IndependentProposalKernel* initial_proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                                     evaluate_log_proposal_in);
   
  for (auto i=likelihood_estimators.begin();
            i!=likelihood_estimators.end();
            ++i)
  {
      (*i)->setup(initial_proposal->independent_simulate(*this->rng));
  }
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(initial_proposal,
                                                                                this->model_and_algorithm.likelihood_estimators);

  if (parallel_in==TRUE)
  {
      //this->the_worker = new RcppParallelSMCWorker(this,
                                                //this->model_and_algorithm.particle_simulator,
                                                //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this,
                                               this->model_and_algorithm.particle_simulator,
                                               this->model_and_algorithm.likelihood_estimators);
  }
  
  //std::string variable_in = "power";
  this->model_and_algorithm.smc_criterion = NULL;//new PositiveSMCCriterion();
  this->model_and_algorithm.smc_termination = NULL;
  //this->sequencer = Sequencer(this->the_worker,
  //                            temperatures_in,
  //                            variable_in,
  //                            this->model_and_algorithm.smc_criterion,
  //                            this->model_and_algorithm.smc_termination);
  
  this->proposal_kernel = proposal_kernel_in;
  */
}

//Copy constructor for the ParticleFilter class.
ParticleFilter::ParticleFilter(const ParticleFilter &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the ParticleFilter class.
ParticleFilter::~ParticleFilter(void)
{
  if (this->proposal_kernel!=NULL)
    delete this->proposal_kernel;
}

void ParticleFilter::operator=(const ParticleFilter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal_kernel!=NULL)
    delete this->proposal_kernel;

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
  if (another.proposal_kernel!=NULL)
    this->proposal_kernel = another.proposal_kernel;
  else
    this->proposal_kernel = another.proposal_kernel;
  
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
  this->dynamic_proposal_is_evaluated = another.dynamic_proposal_is_evaluated;
}

SMCOutput* ParticleFilter::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

/*
SMCOutput* ParticleFilter::specific_run(const std::string &directory_name)
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}
*/

SMCOutput* ParticleFilter::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->factors->set_data(this->current_index);
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
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state);
    //this->sequencer.find_next_target_bisection(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->move(next_particles,
                           current_particles);
    this->evaluate_smcfixed_part_smc(current_state);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::evaluate_smc(SMCOutput* current_state)
{
  if (!this->last_index_is_fixed)
  {
    Rcpp::stop("ParticleFilter::evaluate - need to read in a parameter to determine last measurement index.");
  }
  
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void ParticleFilter::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  VectorSingleIndex index(this->current_index);
  
  this->the_worker->smcfixed_weight(&index,
                                    current_state->back());
  
  /*
  if (this->sequencer_parameters!=NULL)
  {
    this->the_worker->smcfixed_weight(&index,
                                      current_state->back(),
                                      *this->sequencer_parameters);
  }
  else
  {
    this->the_worker->smcfixed_weight(&index,
                                      current_state->back());
  }
  */
  //current_state->initialise_next_step();
}

MoveOutput* ParticleFilter::move(RandomNumberGenerator &rng,
                                 Particle &particle)
{
  Particle moved_particle;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    if (i==0)
    {
      //VectorSingleIndex index(0);
      moved_particle = this->proposal_kernel->move(rng,
                                                   particle);
    }
    else
    {
      //VectorSingleIndex index(this->current_index);
      moved_particle = this->proposal_kernel->move(rng,
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
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->pf_initial_weight(current_state->back());
    }
    else
    {
      VectorSingleIndex index(this->current_index);
      
      if (this->dynamic_proposal_is_evaluated)
      {
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    this->proposal_kernel);
      }
      else
      {
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    NULL);
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //current_state->back().set_previous_target_evaluated_to_target_evaluated();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->factors->set_data(this->current_index);
    
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

/*
SMCOutput* ParticleFilter::specific_run(const std::string &directory_name,
                                        const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}
*/

SMCOutput* ParticleFilter::specific_initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  this->first_index = 0;
  this->current_index = this->first_index;
  this->factors->set_data(this->current_index);
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
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    //this->the_worker->move(next_particles,
    //                       current_particles,
    //                       conditioned_on_parameters);
    this->the_worker->move(next_particles,
                           current_particles);
    this->evaluate_smcfixed_part_smc(current_state,
                                     conditioned_on_parameters);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }

}

void ParticleFilter::evaluate_smc(SMCOutput* current_state,
                                  const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void ParticleFilter::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                const Parameters &conditioned_on_parameters)
{
  VectorSingleIndex index(this->current_index);
  this->the_worker->smcfixed_weight(&index,
                                    current_state->back());
  
  /*
  if (this->sequencer_parameters!=NULL)
  {
    Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
    this->the_worker->smcfixed_weight(&index,
                                      current_state->back(),
                                      all_parameters);
  }
  else
  {
    this->the_worker->smcfixed_weight(&index,
                                      current_state->back(),
                                      conditioned_on_parameters);
  }
  */
  //current_state->initialise_next_step();
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
   
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->pf_initial_weight(current_state->back());
    }
    else
    {
      VectorSingleIndex index(this->current_index);
      
      if (this->dynamic_proposal_is_evaluated)
      {
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    this->proposal_kernel);
        
        /*
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    this->proposal_kernel,
                                    conditioned_on_parameters);
        */
      }
      else
      {
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    NULL);
        /*
        this->the_worker->pf_weight(&index,
                                    current_state->back(),
                                    *(current_state->end()-2),
                                    NULL,
                                    conditioned_on_parameters);
        */
      }
    }

    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //current_state->back().set_previous_target_evaluated_to_target_evaluated();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->factors->set_data(this->current_index);
    
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
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.subsample_find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    /*
     this->the_worker->subsample_move(next_particles,
     current_particles,
     conditioned_on_parameters);
     */
    this->subsample_evaluate_smcfixed_part_smc(current_state,
                                               conditioned_on_parameters);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::subsample_evaluate_smc(SMCOutput* current_state)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void ParticleFilter::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  VectorSingleIndex index(this->current_index);
  /*
   if (this->sequencer_parameters!=NULL)
   {
   Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
   this->the_worker->subsample_smcfixed_weight(&index,
   current_state->back(),
   all_parameters);
   }
   else
   {
   this->the_worker->subsample_smcfixed_weight(&index,
   current_state->back(),
   conditioned_on_parameters);
   }
   */
  this->the_worker->subsample_smcfixed_weight(&index,
                                              current_state->back());
  //current_state->initialise_next_step();
}

void ParticleFilter::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->subsample_pf_initial_weight(current_state->back());
      /*
       this->the_worker->subsample_pf_initial_weight(current_state->back(),
       conditioned_on_parameters);
       */
      
    }
    else
    {
      VectorSingleIndex index(this->current_index);
      if (this->dynamic_proposal_is_evaluated)
      {
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              this->proposal_kernel);
        /*
         this->the_worker->subsample_pf_weight(&index,
         current_state->back(),
         *(current_state->end()-2),
         this->proposal_kernel,
         conditioned_on_parameters);
         */
      }
      else
      {
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              NULL);
        /*
         this->the_worker->subsample_pf_weight(&index,
         current_state->back(),
         *(current_state->end()-2),
         NULL,
         conditioned_on_parameters);
         */
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->factors->set_data(this->current_index);
    
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
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    //current_state->resample();
    
    // Different at first step?
    //this->sequencer.find_desired_criterion(current_state,
    //                                       conditioned_on_parameters);
    //this->sequencer.subsample_find_next_target_bisection(current_state, conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    /*
    this->the_worker->subsample_move(next_particles,
                                     current_particles,
                                     conditioned_on_parameters);
    */
    this->subsample_evaluate_smcfixed_part_smc(current_state);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void ParticleFilter::subsample_evaluate_smc(SMCOutput* current_state,
                                            const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void ParticleFilter::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                          const Parameters &conditioned_on_parameters)
{
  VectorSingleIndex index(this->current_index);
  /*
  if (this->sequencer_parameters!=NULL)
  {
    Parameters all_parameters = conditioned_on_parameters.merge(*this->sequencer_parameters);
    this->the_worker->subsample_smcfixed_weight(&index,
                                                current_state->back(),
                                                all_parameters);
  }
  else
  {
    this->the_worker->subsample_smcfixed_weight(&index,
                                                current_state->back(),
                                                conditioned_on_parameters);
  }
  */
  this->the_worker->subsample_smcfixed_weight(&index,
                                              current_state->back());
  //current_state->initialise_next_step();
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
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->subsample_pf_initial_weight(current_state->back());
      /*
      this->the_worker->subsample_pf_initial_weight(current_state->back(),
                                                    conditioned_on_parameters);
      */
      
    }
    else
    {
      VectorSingleIndex index(this->current_index);
      if (this->dynamic_proposal_is_evaluated)
      {
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              this->proposal_kernel);
        /*
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              this->proposal_kernel,
                                              conditioned_on_parameters);
        */
      }
      else
      {
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              NULL);
        /*
        this->the_worker->subsample_pf_weight(&index,
                                              current_state->back(),
                                              *(current_state->end()-2),
                                              NULL,
                                              conditioned_on_parameters);
        */
      }
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    
    //if (this->sequencer_parameters!=NULL)
    //  current_state->back().schedule_parameters = *this->sequencer_parameters;
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    //current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->check_termination())
    {
      //terminate = TRUE;
      break;
    }
    
    this->current_index = this->current_index+1;
    this->factors->set_data(this->current_index);
    
    this->subsample_simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

void ParticleFilter::weight_for_adapting_sequence(const Index* index,
                                                  Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_evaluate_target(index,
                                                               current_particles);
}

/*
void ParticleFilter::weight_for_adapting_sequence(const Index* index,
                                                  Particles &current_particles,
                                                  const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_evaluate_target(index,
                                                               current_particles,
                                                               conditioned_on_parameters);
}
*/

MoveOutput* ParticleFilter::subsample_move(RandomNumberGenerator &rng,
                                           Particle &particle)
{
  Particle moved_particle;
  for (size_t i=0; i<this->predictions_per_update; ++i)
  {
    if (i==0)
    {
      //VectorSingleIndex index(0);
      moved_particle = this->proposal_kernel->subsample_move(rng,
                                                             particle);
    }
    else
    {
      //VectorSingleIndex index(this->current_index);
      moved_particle = this->proposal_kernel->subsample_move(rng,
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

/*
void ParticleFilter::subsample_weight_for_adapting_sequence(const Index* index,
                                                            Particles &current_particles,
                                                            const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_evaluate_target(index,
                                                                         current_particles,
                                                                         conditioned_on_parameters);
}
*/

bool ParticleFilter::check_termination() const
{
  return(this->current_index==this->last_index);
}

// void ParticleFilter::smc_step(void)
// {
//
// }
//
// void ParticleFilter::weight_update(void)
// {
//
// }

//void ParticleFilter::smc_update(SMCOutput* current_state)
//{
  /*
   unsigned int number_of_points = algorithm["number_of_points"];

   List observed_data = model["observed_data"];

   // Do the initial importance sampling step.

   // Do the simulation.
   SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
   SimulateDistributionPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

   std::vector<List> proposed_points;
   proposed_points.reserve(number_of_points);
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   proposed_points.push_back(simulate_proposal());
   }

   LikelihoodEstimator* likelihood_estimator = make_likelihood_estimator(model, algorithm);

   std::vector<List> proposed_auxiliary_variables;
   proposed_auxiliary_variables.reserve(number_of_points);

   for (std::vector<List>::const_iterator i=proposed_points.begin(); i!=proposed_points.end(); ++i)
   {
   proposed_auxiliary_variables.push_back(likelihood_estimator->simulate_auxiliary_variables(*i));
   }

   likelihood_estimator->is_setup_likelihood_estimator(proposed_points,
   proposed_auxiliary_variables);

   arma::colvec log_weights(number_of_points);
   bool prior_is_proposal = algorithm["prior_is_proposal"];
   if (prior_is_proposal==TRUE)
   {
   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]);
   }
   }
   else
   {
   SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
   EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

   SEXP evaluate_log_proposal_SEXP = algorithm["evaluate_log_proposal"];
   EvaluateLogDistributionPtr evaluate_log_proposal = load_evaluate_log_distribution(evaluate_log_proposal_SEXP);

   for (unsigned int i=0; i<number_of_points; ++i)
   {
   log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]) + evaluate_log_prior(proposed_points[i]) - evaluate_log_proposal(proposed_points[i]);
   }
   }

   if (likelihood_estimator != NULL)
   delete likelihood_estimator;

   return List::create(Named("proposed_points") = proposed_points,
   Named("proposed_auxiliary_variables") = wrap(proposed_auxiliary_variables),
   Named("log_weights") = log_weights,
   Named("log_normalising_constant") = log_sum_exp(log_weights));
   */

//}
