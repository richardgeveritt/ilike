#include <numeric>
#include "smc_marginal.h"
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
#include "vector_factors.h"
#include "custom_distribution_proposal_kernel.h"
#include "positive_smc_criterion.h"
#include "vector_single_index.h"

SMCMarginal::SMCMarginal()
   :SMC()
{
  this->smc_criterion = NULL;
  this->index = NULL;
}

SMCMarginal::SMCMarginal(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         ProposalKernel* proposal_kernel_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateIndependentProposalPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         const std::vector<double> &temperatures_in,
                         bool parallel_in,
                         size_t grain_size_in)
  :SMC(rng_in, seed_in, data_in, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, double(number_of_particles_in), true, false, true)
{
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                               seed_in,
                               data_in,
                               evaluate_log_prior_in,
                               evaluate_log_likelihood_in,
                               TRUE);
  
  EvaluateLogDistributionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                  seed_in,
                                  data_in,
                                  prior_times_likelihood,
                                  power,
                                  TRUE));
  indices.push_back(0);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                               seed_in,
                               data_in,
                               proposal,
                               TRUE);
  
  EvaluateLogDistributionPtr second_power = annealing_one_minus_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  TRUE));
  indices.push_back(1);
  this->index = new VectorSingleIndex(indices);
                                         
  this->factors = new VectorFactors(likelihood_estimators);
  
  

  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);

  if (parallel_in==TRUE)
  {
      //this->the_worker = new RcppParallelSMCWorker(this,
                                                //this->model_and_algorithm.particle_simulator,
                                                //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::string variable_in = "power";
  this->smc_criterion = new PositiveSMCCriterion();
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              temperatures_in,
                              variable_in,
                              this->smc_criterion,
                              this->smc_termination);
  
  this->proposal_kernel = proposal_kernel_in;
}

//Copy constructor for the SMCMarginal class.
SMCMarginal::SMCMarginal(const SMCMarginal &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCMarginal class.
SMCMarginal::~SMCMarginal(void)
{
  if (this->proposal_kernel!=NULL)
    delete this->proposal_kernel;
  
  if (this->smc_criterion!=NULL)
    delete this->smc_criterion;
  
  if (this->index!=NULL)
    delete this->index;
}

void SMCMarginal::operator=(const SMCMarginal &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->proposal_kernel!=NULL)
    delete this->proposal_kernel;
  
  if (this->smc_criterion!=NULL)
    delete this->smc_criterion;
  
  if (this->index!=NULL)
    delete this->index;

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCMarginal::smc_duplicate(void) const
{
  return( new SMCMarginal(*this));
}

LikelihoodEstimator* SMCMarginal::duplicate(void) const
{
  return( new SMCMarginal(*this));
}

void SMCMarginal::make_copy(const SMCMarginal &another)
{
  if (another.proposal_kernel!=NULL)
    this->proposal_kernel = another.proposal_kernel->proposal_kernel_duplicate();
  else
    this->proposal_kernel = NULL;
  
  if (another.smc_criterion!=NULL)
    this->smc_criterion = another.smc_criterion->duplicate();
  else
    this->smc_criterion = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

SMCOutput* SMCMarginal::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* SMCMarginal::initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void SMCMarginal::simulate_smc(SMCOutput* current_state)
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
    
    current_state->normalise_weights();
    current_state->resample();
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state);
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->move(next_particles,
                           current_particles);
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMarginal::evaluate_smc(SMCOutput* current_state)
{
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMarginal::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
}

void SMCMarginal::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    std::runtime_error("SMCMarginal::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->weight(this->index,
                               current_state->back());
      
    }
    else
    {
      this->the_worker->marginal_weight(this->index,
                                        current_state->back(),
                                        *(current_state->end()-2),
                                        proposal_kernel);
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->simulate_smc(current_state);
    
  }
}


// write different types of weight in worker
// weight update is specific to particular SMC
// choosing desired criterion - use ratio to choose desired ess/cess
// finding next target -> weight here appears to be related to type of SMC

MoveOutput* SMCMarginal::move(RandomNumberGenerator &rng,
                              Particle &particle)
{
  return new SinglePointMoveOutput(this->proposal_kernel->move(rng,
                                                               particle));
}

void SMCMarginal::weight_for_adapting_sequence(Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_evaluate_target(this->index,
                                                               current_particles);
}


SMCOutput* SMCMarginal::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* SMCMarginal::initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void SMCMarginal::simulate_smc(SMCOutput* current_state,
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
    
    current_state->normalise_weights();
    current_state->resample();
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state,
                                           conditioned_on_parameters);
    this->sequencer.find_next_target_bisection(current_state,
                                               this->index,
                                               conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->move(next_particles,
                           current_particles,
                           conditioned_on_parameters);
    // involves complete evaluation of weights using current adaptive param
  }

}

void SMCMarginal::evaluate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCMarginal::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                          const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back(),
                                    conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void SMCMarginal::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                               const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
   
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->weight(this->index,
                               current_state->back(),
                               conditioned_on_parameters);
      
    }
    else
    {
      this->the_worker->marginal_weight(this->index,
                                        current_state->back(),
                                        *(current_state->end()-2),
                                        proposal_kernel,
                                        conditioned_on_parameters);
    }

    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

void SMCMarginal::subsample_simulate_smc(SMCOutput* current_state,
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
    
    current_state->normalise_weights();
    current_state->resample();
    
    // Different at first step?
    this->sequencer.find_desired_criterion(current_state,
                                           conditioned_on_parameters);
    this->sequencer.subsample_find_next_target_bisection(current_state,
                                                         this->index,
                                                         conditioned_on_parameters);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->subsample_move(next_particles,
                           current_particles,
                           conditioned_on_parameters);
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMarginal::subsample_evaluate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                             conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                               conditioned_on_parameters);
}

void SMCMarginal::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcfixed_weight(this->index,
                                              current_state->back(),
                                              conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void SMCMarginal::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                         const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // If first step, then weight then check termination. What to do about find desired criterion and find next target? See below - make it different in first iteration?
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    if (current_state->number_of_smc_iterations()==1)
    {
      this->the_worker->subsample_weight(this->index,
                                         current_state->back(),
                                         conditioned_on_parameters);
      
    }
    else
    {
      this->the_worker->subsample_marginal_weight(this->index,
                                                  current_state->back(),
                                                  *(current_state->end()-2),
                                                  proposal_kernel,
                                                  conditioned_on_parameters);
    }
    
    current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->subsample_simulate_smc(current_state, conditioned_on_parameters);
    
  }
}

// write different types of weight in worker
// weight update is specific to particular SMC
// choosing desired criterion - use ratio to choose desired ess/cess
// finding next target -> weight here appears to be related to type of SMC

MoveOutput* SMCMarginal::move(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters)
{
  return new SinglePointMoveOutput(this->proposal_kernel->move(rng,
                                                               particle,
                                                               conditioned_on_parameters));
}

void SMCMarginal::weight_for_adapting_sequence(Particles &current_particles,
                                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_evaluate_target(this->index,
                                                               current_particles,
                                                               conditioned_on_parameters);
}

MoveOutput* SMCMarginal::subsample_move(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters)
{
  return new SinglePointMoveOutput(this->proposal_kernel->subsample_move(rng,
                                                                         particle,
                                                                         conditioned_on_parameters));
}

void SMCMarginal::subsample_weight_for_adapting_sequence(Particles &current_particles,
                                                         const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_evaluate_target(this->index,
                                                                         current_particles,
                                                                         conditioned_on_parameters);
}

// void SMCMarginal::smc_step(void)
// {
//
// }
//
// void SMCMarginal::weight_update(void)
// {
//
// }

//void SMCMarginal::smc_update(SMCOutput* current_state)
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
