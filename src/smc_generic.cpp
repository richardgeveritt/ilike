#include <numeric>
#include "smc_generic.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "positive_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "vector_factors.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_index.h"

namespace ilike
{
SMCGeneric::SMCGeneric()
:SMC()
{
  this->index.reset();
}

SMCGeneric::SMCGeneric(RandomNumberGenerator* rng_in,
                       size_t* seed_in,
                       Data* data_in,
                       size_t number_of_particles_in,
                       size_t lag_in,
                       size_t lag_proposed_in,
                       double resampling_desired_ess_in,
                       const std::vector<double> &temperatures_in,
                       ProposalKernel* proposal_kernel_in,
                       ProposalKernel* L_kernel_in,
                       EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                       EvaluateLogDistributionPtr evaluate_log_prior_in,
                       SimulateDistributionPtr simulate_proposal_in,
                       EvaluateLogDistributionPtr evaluate_log_proposal_in,
                       bool transform_proposed_particles,
                       bool parallel_in,
                       size_t grain_size_in,
                       const std::string &results_name_in)
:SMC(rng_in, seed_in, data_in, Parameters(), number_of_particles_in, lag_in, lag_proposed_in, proposal_kernel_in->get_proposals(), resampling_desired_ess_in, true, false, true, transform_proposed_particles, results_name_in)
{
  proposal_kernel_in->set_proposal_parameters(&this->algorithm_parameters);
  L_kernel_in->set_proposal_parameters(&this->algorithm_parameters);
  
  std::string variable_in = "power";
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  proposal->set_proposal_parameters(&this->algorithm_parameters);

  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  TRUE);
  
  PowerFunctionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  prior_times_likelihood,
                                                                  power,
                                                                  variable_in,
                                                                  true));
  indices.push_back(0);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   evaluate_log_proposal_in,
                                                                                   TRUE);
  
  PowerFunctionPtr second_power = annealing_one_minus_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  variable_in,
                                                                  true));
  indices.push_back(1);
  this->index.reset(new VectorIndex(indices));

  this->factors.reset(new VectorFactors(likelihood_estimators));
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator.reset(new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators));
  
  if (parallel_in==TRUE)
  {
  }
  else
  {
    this->the_worker.reset(new SequentialSMCWorker(this));
  }
  
  SMCCriterion* smc_criterion = new PositiveSMCCriterion();
  this->sequencer = Sequencer(this->the_worker.get(),
                              temperatures_in,
                              variable_in,
                              25,
                              smc_criterion);
  
  this->proposal_kernel.reset(proposal_kernel_in);
  this->L_kernel.reset(L_kernel_in);
}

//Copy constructor for the SMCGeneric class.
SMCGeneric::SMCGeneric(const SMCGeneric &another)
:SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCGeneric class.
SMCGeneric::~SMCGeneric()
{
  this->proposal_kernel.reset();
  
  this->L_kernel.reset();
  
  this->index.reset();
}

void SMCGeneric::operator=(const SMCGeneric &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal_kernel.reset();
  
  this->L_kernel.reset();
  
  this->index.reset();
  
  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCGeneric::smc_duplicate() const
{
  return( new SMCGeneric(*this));
}

LikelihoodEstimator* SMCGeneric::duplicate() const
{
  return( new SMCGeneric(*this));
}

void SMCGeneric::make_copy(const SMCGeneric &another)
{
  if (another.proposal_kernel!=NULL)
    this->proposal_kernel.reset(another.proposal_kernel->proposal_kernel_duplicate());
  else
    this->proposal_kernel.reset();
  
  if (another.L_kernel!=NULL)
    this->L_kernel.reset(another.L_kernel->proposal_kernel_duplicate());
  else
    this->L_kernel.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

SMCOutput* SMCGeneric::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* SMCGeneric::specific_initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

void SMCGeneric::simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state);
    
    // at some point this one might need to be changed to bisection on the ESS
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index.get());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
  }
  else
  {
    // update proposals
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state);
    this->sequencer.find_next_target_quantile(current_state,
                                              index.get());
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->move(next_particles,
                           current_particles);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    this->evaluate_smcfixed_part_smc(current_state);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCGeneric::evaluate_smc(SMCOutput* current_state)
{
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCGeneric::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
}

void SMCGeneric::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    Rcpp::stop("SMCGeneric::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    if (current_state->number_of_smc_iterations()==1)
    {
      
      this->the_worker->weight(this->index.get(),
                               current_state->back());
    }
    else
    {
      this->L_kernel->smc_adapt(current_state);

      this->the_worker->generic_weight(this->index.get(),
                                       current_state->back(),
                                       *(current_state->end()-2),
                                       this->proposal_kernel.get(),
                                       this->L_kernel.get());
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);

    current_state->set_time_and_reset_start();
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state);
    
  }
}

MoveOutput* SMCGeneric::move(RandomNumberGenerator &rng,
                             const Particle &particle) const
{
  return new SinglePointMoveOutput(this->proposal_kernel->move(rng,
                                                               particle));
}

void SMCGeneric::weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index.get(),
                                                      current_particles);
}

SMCOutput* SMCGeneric::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_and_resample_weights();
  return simulation;
}

SMCOutput* SMCGeneric::specific_initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed, this->results_name);
  return output;
}

void SMCGeneric::simulate_smc(SMCOutput* current_state,
                              const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state);
    
    // at some point this one might need to be changed to bisection on the ESS
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index.get());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
  }
  else
  {
    // update proposals
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state);
    this->sequencer.find_next_target_quantile(current_state,
                                              this->index.get());
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    
    this->the_worker->move(next_particles,
                           current_particles);

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    this->evaluate_smcfixed_part_smc(current_state,
                                     conditioned_on_parameters);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCGeneric::evaluate_smc(SMCOutput* current_state,
                              const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCGeneric::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                            const Parameters &conditioned_on_parameters)
{
  
  this->the_worker->smcfixed_weight(this->index.get(),
                                    current_state->back());
}

void SMCGeneric::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
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
      this->the_worker->weight(this->index.get(),
                               current_state->back());
    }
    else
    {
      this->L_kernel->smc_adapt(current_state);
      
      this->the_worker->generic_weight(this->index.get(),
                                       current_state->back(),
                                       *(current_state->end()-2),
                                       this->proposal_kernel.get(),
                                       this->L_kernel.get());
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    current_state->set_llhd(current_state->log_likelihood);
    
    current_state->set_time_and_reset_start();

    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

MoveOutput* SMCGeneric::subsample_move(RandomNumberGenerator &rng,
                                       const Particle &particle) const
{
  return new SinglePointMoveOutput(this->proposal_kernel->subsample_move(rng,
                                                                         particle));
}

void SMCGeneric::subsample_weight_for_adapting_sequence(const Index* index,
                                                        Particles &current_particles)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index.get(),
                                                                current_particles);
}

void SMCGeneric::subsample_simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
    
    // Different at first step?
    this->sequencer.subsample_find_desired_criterion(current_state);
    
    // at some point this one might need to be changed to bisection on the ESS
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index.get());

    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
  }
  else
  {
    // update proposals
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    this->sequencer.subsample_find_desired_criterion(current_state);
    this->sequencer.subsample_find_next_target_quantile(current_state,
                                                        this->index.get());
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    this->subsample_evaluate_smcfixed_part_smc(current_state);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCGeneric::subsample_simulate_smc(SMCOutput* current_state,
                                        const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
    
    // Different at first step?
    this->sequencer.subsample_find_desired_criterion(current_state);
    
    // at some point this one might need to be changed to bisection on the ESS
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index.get());
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
  }
  else
  {
    // update proposals
    this->proposal_kernel->smc_adapt(current_state);
    
    current_state->normalise_and_resample_weights();
    
    // Different at first step?
    this->sequencer.subsample_find_desired_criterion(current_state);
    this->sequencer.subsample_find_next_target_quantile(current_state,
                                                        this->index.get());
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles(current_particles);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();
    
    this->subsample_evaluate_smcfixed_part_smc(current_state,
                                               conditioned_on_parameters);
    
    current_state->increment_smc_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCGeneric::subsample_evaluate_smc(SMCOutput* current_state)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCGeneric::subsample_evaluate_smc(SMCOutput* current_state,
                                        const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void SMCGeneric::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void SMCGeneric::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                      const Parameters &conditioned_on_parameters)
{
  
  this->the_worker->subsample_smcfixed_weight(this->index.get(),
                                              current_state->back());
}

void SMCGeneric::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->subsample_weight(this->index.get(),
                                         current_state->back());
    }
    else
    {
      this->L_kernel->smc_adapt(current_state);
      
      this->the_worker->subsample_generic_weight(this->index.get(),
                                                 current_state->back(),
                                                 *(current_state->end()-2),
                                                 this->proposal_kernel.get(),
                                                 this->L_kernel.get());
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state);
    
  }
}

void SMCGeneric::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
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
      this->the_worker->subsample_weight(this->index.get(),
                                         current_state->back());
    }
    else
    {
      this->L_kernel->smc_adapt(current_state);
      
      this->the_worker->subsample_generic_weight(this->index.get(),
                                                 current_state->back(),
                                                 *(current_state->end()-2),
                                                 this->proposal_kernel.get(),
                                                 this->L_kernel.get());
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();

    current_state->set_llhd(current_state->subsample_log_likelihood);
    current_state->set_time_and_reset_start();
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      current_state->terminate();
      break;
    }
    
    this->subsample_simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

}
