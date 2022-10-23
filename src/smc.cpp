#include "smc.h"

//#include "function_pointers.h"
#include "utils.h"
//#include "likelihood_estimator.h"
//#include "likelihood_maker.h"
#include "smc_worker.h"
#include "sequential_smc_worker.h"
#include "smc_output.h"
#include "smc_criterion.h"
#include "ess_smc_criterion.h"
#include "move_output.h"
#include "particle_simulator.h"
#include "smc_termination.h"
#include "ess_smc_criterion.h"
#include "factor_variables.h"

SMC::SMC()
  :LikelihoodEstimator()
{
  this->particle_simulator = NULL;
  this->smc_termination = NULL;
  //this->output = NULL;
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         Data* data_in,
         size_t number_of_particles_in,
         size_t lag_in,
         size_t lag_proposed_in,
         double resampling_desired_ess_in,
         bool proposal_is_evaluated_in,
         //EvaluateLogDistributionPtr evaluate_log_proposal_in,
         bool smcfixed_flag_in,
         bool sequencer_limit_is_fixed_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in)
{
  //this->output = new SMCOutput(lag_in,
  //                             lag_proposed_in);
  this->number_of_particles = number_of_particles_in;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->smcfixed_flag = smcfixed_flag_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  //this->evaluate_log_proposal = evaluate_log_proposal_in;
  this->proposal_is_evaluated = proposal_is_evaluated_in;
  this->resampling_criterion = new ESSSMCCriterion(resampling_desired_ess_in);
  this->proposed_particles_inputted = false;
  this->initial_particles = Particles();
  this->particle_simulator = NULL;
  this->smc_termination = NULL;
  // Set up worker?
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         Data* data_in,
         size_t lag_in,
         size_t lag_proposed_in,
         double resampling_desired_ess_in,
         const std::vector<Parameters> &initial_values_in,
         const arma::colvec &log_probabilities_of_initial_values_in,
         bool smcfixed_flag_in,
         bool sequencer_limit_is_fixed_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->number_of_particles = initial_values_in.size();
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  this->smcfixed_flag = smcfixed_flag_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  this->resampling_criterion = new ESSSMCCriterion(resampling_desired_ess_in);
  //this->evaluate_log_proposal = NULL;
  this->proposal_is_evaluated = false;
  this->proposed_particles_inputted = true;
  this->initial_particles = Particles(initial_values_in,
                                      log_probabilities_of_initial_values_in);
  this->particle_simulator = NULL;
  this->smc_termination = NULL;
}

SMC::SMC(const SMC &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

SMC::~SMC(void)
{
  if (this->the_worker!=NULL)
    delete this->the_worker;
  
  if (this->resampling_criterion!=NULL)
    delete this->resampling_criterion;
  
  if (this->particle_simulator!=NULL)
    delete this->particle_simulator;
  
  if (this->smc_termination!=NULL)
    delete this->smc_termination;
  
  //if (this->output!=NULL)
  //  delete this->output;
}

void SMC::operator=(const SMC &another)
{
  if(this == &another)
    return;
  
  if (this->the_worker!=NULL)
    delete this->the_worker;
  
  if (this->resampling_criterion!=NULL)
    delete this->resampling_criterion;
  
  if (this->particle_simulator!=NULL)
    delete this->particle_simulator;
  
  if (this->smc_termination!=NULL)
    delete this->smc_termination;

  this->make_copy(another);
}

void SMC::make_copy(const SMC &another)
{
  if (another.the_worker!=NULL)
    this->the_worker = another.the_worker->duplicate();
  else
    this->the_worker = NULL;
  
  if (another.resampling_criterion!=NULL)
    this->resampling_criterion = another.resampling_criterion->duplicate();
  else
    this->resampling_criterion = NULL;
  
  if (another.particle_simulator!=NULL)
    this->particle_simulator = another.particle_simulator->duplicate();
  else
    this->particle_simulator = NULL;
  
  if (another.smc_termination!=NULL)
    this->smc_termination = another.smc_termination->duplicate();
  else
    this->smc_termination = NULL;

  //if (this->output!=NULL)
  //  this->output = another.output->smc_duplicate();
  
  this->sequencer = another.sequencer;
  this->smcfixed_flag = another.smcfixed_flag;
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  this->number_of_particles = another.number_of_particles;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  //this->evaluate_log_proposal = another.evaluate_log_proposal;
  this->proposal_is_evaluated = another.proposal_is_evaluated;
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_particles = another.initial_particles;
}

//Particles SMC::is_step() const
//{
  // The way in which this is done is determined by what is set in model_and_algorithm.
  //this->model_and_algorithm;

  // One choice will use an RcppParallel worker.
//  return Particles();
//}

// Parameters SMC::single_particle_is_step() const
// {
//   //Parameters result = this->model_and_algorithm->simulate_priors->simulate();
//
//   return Parameters();
//
// }

SMCOutput* SMC::run()
{
  return this->specific_run();
}

SMCOutput* SMC::run(const Parameters &conditioned_on_parameters)
{
  return this->specific_run(conditioned_on_parameters);
}

LikelihoodEstimatorOutput* SMC::initialise()
{
  return this->initialise_smc();
}

void SMC::simulate_proposal(SMCOutput* current_state)
{
  // Simulate from proposal.
  Particles* next_particles = current_state->add_particles();
  
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_particles);
    //Particles particles(this->the_worker->get_particles());
    //current_state->initialise_normalised_log_weights();
    //current_state->initialise_unnormalised_log_incremental_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  }
  else
  {
    *next_particles = this->initial_particles;
  }
  
  //std::cout<<next_particles->back()->back().factor_variables->get_particle()->parameters<<std::endl;
  //std::cout<<current_state->back().back()->back().factor_variables->get_particle()->parameters<<std::endl;
  
  current_state->back().initialise_weights();
}

LikelihoodEstimatorOutput* SMC::initialise(const Parameters &conditioned_on_parameters)
{
  return this->initialise_smc(conditioned_on_parameters);
}

void SMC::simulate_proposal(SMCOutput* current_state,
                            const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  Particles* next_particles = current_state->add_particles();
  
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_particles,
                               conditioned_on_parameters);
  }
  else
  {
    *next_particles = this->initial_particles;
  }
  
  current_state->back().initialise_weights();
}

void SMC::resample(SMCOutput* current_state)
{
  // Check criterion to see if we resample.
  if ((*this->resampling_criterion)(current_state->back()))
  {
    // Sample ancestor variables.
    current_state->back().resample();
  }
  else
  {
    // Set ancester variables to be 1:n.
    std::vector<size_t> ancestor_variables_in;
    ancestor_variables_in.reserve(this->number_of_particles);
    for (size_t i=0; i<this->number_of_particles; ++i)
    {
      ancestor_variables_in.push_back(i);
    }
    current_state->back().ancestor_variables = ancestor_variables_in;
  }
}
