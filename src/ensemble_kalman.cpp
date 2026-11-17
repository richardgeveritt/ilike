#include "ensemble_kalman.h"
#include "ensemble_kalman_output.h"
#include "ensemble_kalman_worker.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "independent_proposal_kernel.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_shifter.h"
#include "ensemble_sequencer.h"
#include "factors.h"
#include "transform.h"

namespace ilike
{
EnsembleKalman::EnsembleKalman()
:LikelihoodEstimator()
{
  this->proposal.reset();
  this->ensemble_factors.reset();
  this->ensemble_shifter.reset();
  this->proposed_particles_inputted = false;
  this->transform = NULL;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  
  this->the_worker.reset();
}

EnsembleKalman::EnsembleKalman(RandomNumberGenerator* rng_in,
                               size_t* seed_in,
                               Data* data_in,
                               size_t number_of_ensemble_members_in,
                               size_t lag_in,
                               EnsembleShifter* shifter_in,
                               std::shared_ptr<Transform> transform_in,
                               bool smcfixed_flag_in,
                               bool sequencer_limit_is_fixed_in,
                               const std::string &results_name_in)
:LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->number_of_ensemble_members = number_of_ensemble_members_in;
  if (this->number_of_ensemble_members<=1)
    Rcpp::stop("EnsembleKalman::EnsembleKalman - number of ensemble members must be more than 1.");
  
  this->lag = lag_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  
  this->proposal.reset();
  this->ensemble_factors.reset();
  this->ensemble_shifter.reset(shifter_in);
  
  this->results_name = results_name_in;
  
  this->proposed_particles_inputted = false;
  
  this->transform = transform_in;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  
  this->the_worker.reset();
}

EnsembleKalman::~EnsembleKalman()
{
  
  this->proposal.reset();
  
  this->ensemble_factors.reset();
  
  this->ensemble_shifter.reset();
  
  this->the_worker.reset();
  
}

//Copy constructor for the EnsembleKalman class.
EnsembleKalman::EnsembleKalman(const EnsembleKalman &another)
:LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void EnsembleKalman::operator=(const EnsembleKalman &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal.reset();
  
  this->ensemble_factors.reset();
  
  this->ensemble_shifter.reset();
  
  this->the_worker.reset();

  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

void EnsembleKalman::make_copy(const EnsembleKalman &another)
{
  this->number_of_ensemble_members = another.number_of_ensemble_members;
  this->lag = another.lag;
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  
  this->reciprocal_schedule_scale = another.reciprocal_schedule_scale;

  if (another.proposal!=NULL)
    this->proposal.reset(another.proposal->independent_proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.ensemble_factors!=NULL)
    this->ensemble_factors.reset(another.ensemble_factors->duplicate());
  else
    this->ensemble_factors.reset();
  
  if (another.ensemble_shifter!=NULL)
    this->ensemble_shifter.reset(another.ensemble_shifter->duplicate());
  else
    this->ensemble_shifter.reset();
  
  if (another.the_worker!=NULL)
    this->the_worker.reset(another.the_worker->duplicate());
  else
    this->the_worker.reset();
  this->the_worker->set_enk(this);
  
  this->transform = another.transform;
  
  this->packing_instructions = another.packing_instructions;
  this->vector_variables = another.vector_variables;
  this->any_variables = another.any_variables;
  this->vector_variable_sizes = another.vector_variable_sizes;
  this->results_name = another.results_name;
  
  this->sequencer = another.sequencer;
  
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_ensemble = another.initial_ensemble;
  
  this->initialised = another.initialised;
  
}

EnsembleKalmanOutput* EnsembleKalman::run()
{
  this->setup();
  return this->specific_run();
}

EnsembleKalmanOutput* EnsembleKalman::run(const Parameters &conditioned_on_parameters)
{
  this->setup(conditioned_on_parameters);
  return this->specific_run(conditioned_on_parameters);
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise()
{
  return this->ensemble_kalman_initialise();
}

EnsembleKalmanOutput* EnsembleKalman::ensemble_kalman_initialise()
{
  if (this->initialised==false)
  {
    this->setup();
    this->initialised = true;
  }
  
  return this->specific_ensemble_kalman_initialise();
}

void EnsembleKalman::setup()
{
  this->setup_variables();
}

void EnsembleKalman::setup(const Parameters &parameters)
{
  this->setup_variables(parameters);
}

void EnsembleKalman::setup_variables()
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
  this->any_variables = dummy_parameters.get_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
  
  if (this->ensemble_factors!=NULL)
    this->ensemble_factors->setup(dummy_parameters);
}

void EnsembleKalman::setup_variables(const Parameters &parameters)
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng,
                                                            parameters);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
  this->any_variables = dummy_parameters.get_any_variables();
  
  this->vector_variable_sizes = dummy_parameters.get_variable_n_elems(this->vector_variables);
  
  dummy_parameters.merge_with_fixed(parameters);
  dummy_parameters.merge_with_fixed(this->sequencer.schedule_parameters);
  if (this->factors!=NULL)
    this->factors->setup(dummy_parameters);
  
  if (this->ensemble_factors!=NULL)
    this->ensemble_factors->setup(dummy_parameters);
}

void EnsembleKalman::set_packing_instructions()
{
  Parameters dummy_parameters;
  if (this->proposed_particles_inputted)
  {
    dummy_parameters = this->initial_ensemble[0];
  }
  else
  {
    dummy_parameters = this->proposal->independent_simulate(*this->rng);
  }
  
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->packing_instructions.set_info(dummy_parameters,this->vector_variables);
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise(const Parameters &conditioned_on_parameters)
{
  return this->ensemble_kalman_initialise(conditioned_on_parameters);
}

EnsembleKalmanOutput* EnsembleKalman::ensemble_kalman_initialise(const Parameters &conditioned_on_parameters)
{
  
  if (this->initialised==false)
  {
    this->setup(conditioned_on_parameters);
    this->initialised = true;
  }
  return this->specific_ensemble_kalman_initialise(conditioned_on_parameters);
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index)
{
  // Simulate from proposal.
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
  
  if (!this->proposed_particles_inputted)
  {
    
    this->the_worker->simulate(next_ensemble,
                               index);
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble,
                         this->ensemble_factors.get(),
                         this->sequencer.schedule_parameters);
  }
  
  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  //Ensemble* next_ensemble = current_state->add_ensemble();
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors.get());
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_ensemble,
                               index,
                               conditioned_on_parameters);
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble,
                         this->ensemble_factors.get(),
                         conditioned_on_parameters,
                         this->sequencer.schedule_parameters);
  }
  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle,
                                              const Parameters &sequencer_parameters) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng),
                        this->ensemble_factors.get(),
                        sequencer_parameters);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng)),
                        this->ensemble_factors.get(),
                        sequencer_parameters);
  }

  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors.get());
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle,
                                              const Parameters &sequencer_parameters,
                                              const Parameters &conditioned_on_parameters) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng,
                                                             conditioned_on_parameters),
                        this->ensemble_factors.get(),
                        conditioned_on_parameters,
                        sequencer_parameters);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng,
                                                                                        conditioned_on_parameters)),
                        this->ensemble_factors.get(),
                        conditioned_on_parameters,
                        sequencer_parameters);
  }

  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors.get(), conditioned_on_parameters);
}

void EnsembleKalman::find_measurement_covariances(EnsembleKalmanOutput* simulation)
{
  simulation->back().find_measurement_covariances();
}

void EnsembleKalman::set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in)
{
  this->reciprocal_schedule_scale = reciprocal_schedule_scale_in;
}

void EnsembleKalman::specific_change_data(Data* new_data)
{
  
}
}
