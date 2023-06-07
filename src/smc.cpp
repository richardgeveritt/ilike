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
#include "transform.h"
#include "factors.h"

SMC::SMC()
  :LikelihoodEstimator()
{
  this->particle_simulator = NULL;
  //this->smc_termination = NULL;
  //this->sequencer_parameters = NULL;
  //this->transform = NULL;
  //this->store_raw = true;
  //this->store_transformed = false;
  this->setup_default_ancestor_variables();
  this->initialised = false;
  this->the_worker = NULL;
  this->resampling_criterion = NULL;
  //this->output = NULL;
}

SMC::SMC(RandomNumberGenerator* rng_in,
         size_t* seed_in,
         Data* data_in,
         const Parameters &algorithm_parameters_in,
         size_t number_of_particles_in,
         size_t lag_in,
         size_t lag_proposed_in,
         double resampling_desired_ess_in,
         bool proposal_is_evaluated_in,
         //EvaluateLogDistributionPtr evaluate_log_proposal_in,
         bool smcfixed_flag_in,
         bool sequencer_limit_is_fixed_in,
         const std::string &results_name_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in, algorithm_parameters_in, smcfixed_flag_in)
{
  //this->output = new SMCOutput(lag_in,
  //                             lag_proposed_in);
  this->number_of_particles = number_of_particles_in;
  this->lag = lag_in;
  this->lag_proposed = lag_proposed_in;
  //this->summary_statistics = summary_statistics_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  //this->evaluate_log_proposal = evaluate_log_proposal_in;
  this->proposal_is_evaluated = proposal_is_evaluated_in;
  this->resampling_criterion = new ESSSMCCriterion(resampling_desired_ess_in);
  this->proposed_particles_inputted = false;
  this->particle_simulator = NULL;
  //this->smc_termination = NULL;
  //this->sequencer_parameters = NULL;
  //this->transform = NULL;
  //this->store_raw = true;
  //this->store_transformed = false;
  this->setup_default_ancestor_variables();
  this->results_name = results_name_in;
  
  this->initialised = false;
  // Set up worker?
}

/*
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
  //this->results_name = results_name_in;
}
*/

SMC::SMC(const SMC &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

SMC::~SMC()
{
  if (this->the_worker!=NULL)
    delete this->the_worker;
  
  if (this->resampling_criterion!=NULL)
    delete this->resampling_criterion;
  
  if (this->particle_simulator!=NULL)
    delete this->particle_simulator;
  
  //if (this->transform!=NULL)
  //  delete this->transform;
  
  //if (this->smc_termination!=NULL)
  //  delete this->smc_termination;
  
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
  
  //if (this->transform!=NULL)
  //  delete this->transform;
  
  //if (this->smc_termination!=NULL)
  //  delete this->smc_termination;

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
  
  //if (another.transform!=NULL)
  //  this->transform = another.transform->duplicate();
  //else
  //  this->transform = NULL;
  
  //if (another.smc_termination!=NULL)
  //  this->smc_termination = another.smc_termination->duplicate();
  //else
  //  this->smc_termination = NULL;

  //if (this->output!=NULL)
  //  this->output = another.output->smc_duplicate();
  
  this->sequencer = another.sequencer;
  /*
  if (another.sequencer_parameters!=NULL)
  {
    this->sequencer_parameters = &this->sequencer.schedule_parameters;
  }
  else
  {
    this->sequencer_parameters = NULL;
  }
  */
  
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  this->number_of_particles = another.number_of_particles;
  this->lag = another.lag;
  this->lag_proposed = another.lag_proposed;
  //this->summary_statistics = another.summary_statistics;
  //this->evaluate_log_proposal = another.evaluate_log_proposal;
  this->proposal_is_evaluated = another.proposal_is_evaluated;
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_particles = another.initial_particles;
  this->vector_variables = another.vector_variables;
  this->any_variables = another.any_variables;
  this->vector_variable_sizes = another.vector_variable_sizes;
  this->default_ancestor_variables = another.default_ancestor_variables;
  //this->store_raw = another.store_raw;
  //this->store_transformed = another.store_transformed;
  //this->output_filename = another.output_filename;
  this->log_probabilities_of_initial_values = another.log_probabilities_of_initial_values;
  this->initialised = another.initialised;
  this->results_name = another.results_name;
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
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  return this->specific_run();
}

/*
SMCOutput* SMC::run(const std::string &directory_name)
{
  return this->specific_run(directory_name);
}
*/

SMCOutput* SMC::run(const Parameters &conditioned_on_parameters)
{
  this->setup(conditioned_on_parameters);
  return this->specific_run(conditioned_on_parameters);
}

/*
SMCOutput* SMC::run(const std::string &directory_name,
                    const Parameters &conditioned_on_parameters)
{
  return this->specific_run(directory_name,
                            conditioned_on_parameters);
}
*/

LikelihoodEstimatorOutput* SMC::initialise()
{
  return this->initialise_smc();
}

SMCOutput* SMC::initialise_smc()
{
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  
  return this->specific_initialise_smc();
}

void SMC::setup()
{
  this->setup_variables();
}

void SMC::setup(const Parameters &parameters)
{
  this->setup_variables(parameters);
}

void SMC::setup_variables()
{
  
  Parameters dummy_parameters;
  
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_particles[0];
    dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  }
  else
  {
    dummy_parameters = std::move(this->particle_simulator->simulate(*this->rng,this->factors,this->sequencer.schedule_parameters).parameters);
  }
  
  this->vector_variables = dummy_parameters.get_nonfixed_vector_variables();
  this->any_variables = dummy_parameters.get_nonfixed_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
}

void SMC::setup_variables(const Parameters &parameters)
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_particles[0];
    dummy_parameters.merge_with_fixed(parameters);
    dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  }
  else
  {
    dummy_parameters = std::move(this->particle_simulator->simulate(*this->rng,
                                                                    this->factors,
                                                                    parameters,
                                                                    this->sequencer.schedule_parameters).parameters);
  }
  this->vector_variables = dummy_parameters.get_nonfixed_vector_variables();
  this->any_variables = dummy_parameters.get_nonfixed_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
}

  /*
  Parameters dummy_parameters;
  if (this->transform==NULL)
  {
    this->vector_variables = dummy_parameters.get_vector_variables();
    this->any_variables = dummy_parameters.get_any_variables();
  }
  else
  {
    if (this->store_transformed==false)
    {
      Rcpp::stop("SMC::setup_variables - need to store transformed auxiliary variables if transform is specified.");
    }
    else
    {
      if (this->proposed_particles_inputted)
      {
        dummy_parameters = this->initial_particles[0];
      }
      else
      {
        Particle new_particle;
        this->particle_simulator->simulate(*this->rng,
                                           &new_particle,
                                           this->factors);
        dummy_parameters = new_particle.parameters;
      }
      
      if (this->store_raw==false)
      {
        dummy_parameters = this->transform->transform(dummy_parameters);
        this->vector_variables = dummy_parameters.get_vector_variables();
        this->any_variables = dummy_parameters.get_any_variables();
      }
      else
      {
        dummy_parameters.add_parameters(this->transform->transform(dummy_parameters));
        this->vector_variables = dummy_parameters.get_vector_variables();
        this->any_variables = dummy_parameters.get_any_variables();
      }
      
    }
  }
  */

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
    //arma::mat tau = this->initial_particles[0]["tau"];
    next_particles->setup(this->initial_particles, this->log_probabilities_of_initial_values, this->factors,this->sequencer.schedule_parameters);
  }

  current_state->back().initialise_weights();
}

LikelihoodEstimatorOutput* SMC::initialise(const Parameters &conditioned_on_parameters)
{
  return this->initialise_smc(conditioned_on_parameters);
}

SMCOutput* SMC::initialise_smc(const Parameters &conditioned_on_parameters)
{
  /*
  Particle new_particle;
  Parameters dummy_parameters = this->particle_simulator->simulate(*this->rng,
                                                                   this->factors,
                                                                   conditioned_on_parameters).parameters;
  
  //Parameters dummy_parameters = this->particle_simulator->simulate(*this->rng,
  //                                                                 this->factors,
  //                                                                 conditioned_on_parameters).parameters;
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->any_variables = dummy_parameters.get_any_variables();
  */
  
  if (this->initialised==false)
  {
    this->setup(conditioned_on_parameters);
    this->initialised = true;
  }
  
  return this->specific_initialise_smc(conditioned_on_parameters);
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
    next_particles->setup(this->initial_particles,
                          this->log_probabilities_of_initial_values,
                          this->factors,
                          conditioned_on_parameters,
                          this->sequencer.schedule_parameters);
  }
  
  current_state->back().initialise_weights();
}

void SMC::resample(SMCOutput* current_state)
{
  auto current_particles_iterator = current_state->end()-1;
  
  // Check criterion to see if we resample.
  if ((*this->resampling_criterion)(*current_particles_iterator)<0.0)
  {
    // Sample ancestor variables.
    current_particles_iterator->resampled_flag = true;
    current_particles_iterator->resample();
  }
  else
  {
    current_particles_iterator->ancestor_variables = this->default_ancestor_variables;
  }
}

void SMC::setup_default_ancestor_variables()
{
  // Set ancestor variables to be 1:n.
  this->default_ancestor_variables.clear();
  this->default_ancestor_variables.reserve(this->number_of_particles);
  for (size_t i=0; i<this->number_of_particles; ++i)
  {
    this->default_ancestor_variables.push_back(i);
  }
}
