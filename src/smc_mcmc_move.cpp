#include <numeric>
#include "smc_mcmc_move.h"
#include "smc_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"
#include "parameter_particle_simulator.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_smc_worker.h"
#include "mcmc.h"
#include "move_output.h"
#include "standard_mcmc_output.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_factors.h"
#include "vector_single_index.h"

SMCMCMCMove::SMCMCMCMove()
   :SMC()
{
  this->mcmc = NULL;
  this->smc_criterion = NULL;
  this->index = NULL;
}

// Multiple MCMC chains.
SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         const std::vector<Parameters> &initial_points_in,
                         const arma::colvec &log_probabilities_of_initial_values_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         bool parallel_in,
                         size_t grain_size_in)
:SMC(rng_in,
     seed_in,
     data_in,
     std::max<size_t>(2,lag_in),
     lag_proposed_in,
     -1.0,
     initial_points_in,
     log_probabilities_of_initial_values_in,
     true,
     true)
{
  this->index = NULL;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  likelihood_estimators.reserve(1);
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               evaluate_log_prior_in,
                                                               evaluate_log_likelihood_in,
                                                               true));
  indices.push_back(0);
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  this->particle_simulator = NULL;
  //this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
  //                                                                              this->model_and_algorithm.likelihood_estimators);
  
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
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  this->smc_criterion = new CESSSMCCriterion(-1.0);
  this->smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              this->smc_criterion,
                              this->smc_termination);
  
  this->mcmc = mcmc_in;
}

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateIndependentProposalPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool parallel_in,
                         size_t grain_size_in)
:SMC(rng_in, seed_in, data_in, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, -1.0, true, true, true)
{
  this->index = NULL;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  likelihood_estimators.reserve(1);
  likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                               seed_in,
                                                               data_in,
                                                               evaluate_log_prior_in,
                                                               evaluate_log_likelihood_in,
                                                               true));
  indices.push_back(0);
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  
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
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  this->smc_criterion = new CESSSMCCriterion(-1.0);
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              this->smc_criterion,
                              this->smc_termination);
  
  this->mcmc = mcmc_in;
}

/*
void SMCMCMCMove::set_multiple_mcmc(RandomNumberGenerator* rng_in,
                                    size_t* seed_in,
                                    const Data* data_in,
                                    MCMC* mcmc_in,
                                    SimulateDistributionPtr simulate_proposal_in,
                                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                                    bool parallel_in,
                                    size_t grain_size_in)
{
  this->model_and_algorithm.likelihood_estimators.resize(0);
  this->model_and_algorithm.likelihood_estimators.reserve(1);
  this->model_and_algorithm.likelihood_estimators.push_back(new ExactLikelihoodEstimator(rng_in,
                                                                                         seed_in,
                                                                                         data_in,
                                                                                         evaluate_log_prior_in,
                                                                                         evaluate_log_likelihood_in,
                                                                                         true));
  
  this->model_and_algorithm.particle_simulator = new ParameterParticleSimulator(simulate_proposal_in,
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
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  this->model_and_algorithm.smc_criterion = new CESSSMCCriterion(-1.0);
  this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              this->model_and_algorithm.smc_criterion,
                              this->model_and_algorithm.smc_termination);
  
  this->mcmc = mcmc_in;
}
*/

SMCMCMCMove::SMCMCMCMove(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         size_t number_of_particles_in,
                         size_t lag_in,
                         size_t lag_proposed_in,
                         MCMC* mcmc_in,
                         double resampling_desired_ess_in,
                         double annealing_desired_cess_in,
                         EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                         EvaluateLogDistributionPtr evaluate_log_prior_in,
                         SimulateIndependentProposalPtr simulate_proposal_in,
                         EvaluateLogDistributionPtr evaluate_log_proposal_in,
                         bool parallel_in,
                         size_t grain_size_in)
  :SMC(rng_in, seed_in, data_in, number_of_particles_in, std::max<size_t>(2,lag_in), lag_proposed_in, resampling_desired_ess_in, true, true, true)
{
  this->index = NULL;
  
  std::vector<LikelihoodEstimator*> likelihood_estimators;
  std::vector<size_t> indices;
  
  // Prior times likelihood.
  ExactLikelihoodEstimator* prior_times_likelihood = new ExactLikelihoodEstimator(rng_in,
                                                                                  seed_in,
                                                                                  data_in,
                                                                                  evaluate_log_prior_in,
                                                                                  evaluate_log_likelihood_in,
                                                                                  true);
  
  EvaluateLogDistributionPtr power = annealing_power;
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                  seed_in,
                                  data_in,
                                  prior_times_likelihood,
                                  power,
                                  true));
  indices.push_back(0);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  IndependentProposalKernel* proposal = new CustomDistributionProposalKernel(simulate_proposal_in,
                                                                             evaluate_log_proposal_in);
  
  // Proposal.
  ExactLikelihoodEstimator* proposal_for_evaluation = new ExactLikelihoodEstimator(rng_in,
                                                                                   seed_in,
                                                                                   data_in,
                                                                                   proposal,
                                                                                   true);
  
  EvaluateLogDistributionPtr second_power = annealing_one_minus_power;
                                                            
  
  likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(rng_in,
                                                                  seed_in,
                                                                  data_in,
                                                                  proposal_for_evaluation,
                                                                  second_power,
                                                                  true));
  
  indices.push_back(1);
  this->index = new VectorSingleIndex(indices);
  
  this->factors = new VectorFactors(likelihood_estimators);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  this->particle_simulator = new ParameterParticleSimulator(proposal,
                                                            likelihood_estimators);

  if (parallel_in==true)
  {
      //this->the_worker = new RcppParallelSMCWorker(this,
                                                //this->model_and_algorithm.particle_simulator,
                                                //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialSMCWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  std::string variable_in = "power";
  this->smc_criterion = new CESSSMCCriterion(annealing_desired_cess_in);
  //this->model_and_algorithm.smc_termination = NULL;
  this->sequencer = Sequencer(this->the_worker,
                              schedule_in,
                              variable_in,
                              this->smc_criterion,
                              this->smc_termination);
  
  this->mcmc = mcmc_in;
}

//Copy constructor for the SMCMCMCMove class.
SMCMCMCMove::SMCMCMCMove(const SMCMCMCMove &another)
  :SMC(another)
{
  this->make_copy(another);
}

//Destructor for the SMCMCMCMove class.
SMCMCMCMove::~SMCMCMCMove()
{
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->smc_criterion!=NULL)
    delete this->smc_criterion;
  
  if (this->index!=NULL)
    delete this->index;
}

void SMCMCMCMove::operator=(const SMCMCMCMove &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->smc_criterion!=NULL)
    delete this->smc_criterion;
  
  if (this->index!=NULL)
    delete this->index;

  SMC::operator=(another);
  this->make_copy(another);
}

SMC* SMCMCMCMove::smc_duplicate(void) const
{
  return( new SMCMCMCMove(*this));
}

LikelihoodEstimator* SMCMCMCMove::duplicate() const
{
  return( new SMCMCMCMove(*this));
}

void SMCMCMCMove::make_copy(const SMCMCMCMove &another)
{
  if (another.mcmc!=NULL)
    this->mcmc = another.mcmc->mcmc_duplicate();
  else
    this->mcmc = NULL;
  
  if (another.smc_criterion!=NULL)
    this->smc_criterion = another.smc_criterion->duplicate();
  else
    this->smc_criterion = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

SMCOutput* SMCMCMCMove::specific_run()
{
  SMCOutput* simulation = this->initialise_smc();
  this->simulate_smc(simulation);
  this->evaluate_smc(simulation);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* SMCMCMCMove::initialise_smc()
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_weights();
    current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->move(next_particles,
                           current_particles);
    
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state)
{
  this->evaluate_smcfixed_part_smc(current_state);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state)
{
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back());
  //current_state->initialise_next_step();
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    std::runtime_error("SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    current_state->update_weights(this->sequencer.find_next_target_bisection(current_state,
                                                                             this->index));
    
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

MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              Particle &particle)
{
  return this->mcmc->run(rng,
                         particle);
}

void SMCMCMCMove::weight_for_adapting_sequence(Particles &current_particles)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles);
}

SMCOutput* SMCMCMCMove::specific_run(const Parameters &conditioned_on_parameters)
{
  SMCOutput* simulation = this->initialise_smc(conditioned_on_parameters);
  this->simulate_smc(simulation, conditioned_on_parameters);
  this->evaluate_smc(simulation, conditioned_on_parameters);
  simulation->normalise_weights();
  return simulation;
}

SMCOutput* SMCMCMCMove::initialise_smc(const Parameters &conditioned_on_parameters)
{
  SMCOutput* output = new SMCOutput(this, this->lag, this->lag_proposed);
  return output;
}

void SMCMCMCMove::simulate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_weights();
    current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->move(next_particles,
                           current_particles,
                           conditioned_on_parameters);
    // involves complete evaluation of weights using current adaptive param
  }

}

void SMCMCMCMove::evaluate_smc(SMCOutput* current_state,
                               const Parameters &conditioned_on_parameters)
{
  this->evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCMCMCMove::evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                             const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcfixed_weight(this->index,
                                    current_state->back(),
                                    conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void SMCMCMCMove::evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                               const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state,
                                           conditioned_on_parameters);
    
    // (involves evaluating adaptive weights, using Sequencer)
    current_state->update_weights(this->sequencer.find_next_target_bisection(current_state,
                                                                             this->index,
                                                                             conditioned_on_parameters));
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->simulate_smc(current_state, conditioned_on_parameters);
    
    current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }
}

void SMCMCMCMove::subsample_simulate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_smc_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state, conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    this->mcmc->smc_adapt(current_state);
    
    current_state->normalise_weights();
    current_state->resample();
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Particles* current_particles = &current_state->back();
    Particles* next_particles = current_state->add_particles();
    this->the_worker->subsample_move(next_particles,
                           current_particles,
                           conditioned_on_parameters);
    // involves complete evaluation of weights using current adaptive param
  }
  
}

void SMCMCMCMove::subsample_evaluate_smc(SMCOutput* current_state,
                                         const Parameters &conditioned_on_parameters)
{
  this->subsample_evaluate_smcfixed_part_smc(current_state,
                                   conditioned_on_parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed_smc(current_state,
                                                     conditioned_on_parameters);
}

void SMCMCMCMove::subsample_evaluate_smcfixed_part_smc(SMCOutput* current_state,
                                                       const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcfixed_weight(this->index,
                                              current_state->back(),
                                              conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void SMCMCMCMove::subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* current_state,
                                                                         const Parameters &conditioned_on_parameters)
{
  // set sequencer to have values from conditioned_on_parameters
  if (!this->sequencer_limit_is_fixed)
    this->sequencer.set_next_with_parameter(conditioned_on_parameters);
  
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->sequencer.find_desired_criterion(current_state,
                                           conditioned_on_parameters);
    
    // (involves evaluating adaptive weights, using Sequencer)
    current_state->update_weights(this->sequencer.subsample_find_next_target_bisection(current_state,
                                                                                       this->index,
                                                                                       conditioned_on_parameters));
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->sequencer.check_termination())
    {
      terminate = TRUE;
      break;
    }
    
    this->subsample_simulate_smc(current_state,
                                 conditioned_on_parameters);
    
    current_state->back().subsample_set_previous_target_evaluated_to_target_evaluated();
  }
}

MoveOutput* SMCMCMCMove::move(RandomNumberGenerator &rng,
                              Particle &particle,
                              const Parameters &conditioned_on_parameters)
{
  return this->mcmc->run(rng,
                         particle,
                         conditioned_on_parameters);
}

void SMCMCMCMove::weight_for_adapting_sequence(Particles &current_particles,
                                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(this->index,
                                                      current_particles,
                                                      conditioned_on_parameters);
}

MoveOutput* SMCMCMCMove::subsample_move(RandomNumberGenerator &rng,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters)
{
  return this->mcmc->subsample_run(rng,
                                   particle,
                                   conditioned_on_parameters);
}

void SMCMCMCMove::subsample_weight_for_adapting_sequence(Particles &current_particles,
                                                         const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(this->index,
                                                                current_particles,
                                                                conditioned_on_parameters);
}

// void SMCMCMCMove::smc_step(void)
// {
//
// }
//
// void SMCMCMCMove::weight_update(void)
// {
//
// }

//void SMCMCMCMove::smc_update(SMCOutput* current_state)
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
