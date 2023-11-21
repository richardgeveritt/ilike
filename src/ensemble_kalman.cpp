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

EnsembleKalman::EnsembleKalman()
  :LikelihoodEstimator()
{
  this->proposal = NULL;
  //this->sequencer_parameters = NULL;
  this->ensemble_factors = NULL;
  this->ensemble_shifter = NULL;
  this->proposed_particles_inputted = false;
  this->transform = NULL;
  //this->inverse_transform = NULL;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  //this->measurement_covariance_estimators.resize(0);
  //this->measurements_names.resize(0);
  //this->states_names.resize(0);
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
  this->lag = lag_in;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  
  this->proposal = NULL;
  this->ensemble_factors = NULL;
  this->ensemble_shifter = shifter_in;
  
  //this->sequencer_parameters = NULL;
  this->results_name = results_name_in;
  
  this->proposed_particles_inputted = false;
  
  this->transform = transform_in;
  //this->inverse_transform = inverse_transform_in;
  
  this->initialised = false;
  this->reciprocal_schedule_scale = 0.0;
  //this->measurement_covariance_estimators.resize(0);
  //this->measurements_names.resize(0);
  //this->states_names.resize(0);
  //this->measurements_start_and_end.resize(0);
  //this->states_start_and_end.resize(0);
}

EnsembleKalman::~EnsembleKalman()
{
  /*
  for (std::vector<MeasurementCovarianceEstimator*>::iterator i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  */
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->ensemble_factors!=NULL)
    delete this->ensemble_factors;
  
  if (this->ensemble_shifter!=NULL)
    delete this->ensemble_shifter;
  
  //if (this->transform!=NULL)
  //  delete this->transform;
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
  
  if (this->proposal!=NULL)
    delete this->proposal;
  
  if (this->ensemble_factors!=NULL)
    delete this->ensemble_factors;
  
  if (this->ensemble_shifter!=NULL)
    delete this->ensemble_shifter;
  
  //if (this->transform!=NULL)
  //  delete this->transform;
  
  /*
  for (std::vector<MeasurementCovarianceEstimator*>::iterator i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->measurement_covariance_estimators.clear();
  */

  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

void EnsembleKalman::make_copy(const EnsembleKalman &another)
{
  this->number_of_ensemble_members = another.number_of_ensemble_members;
  this->lag = another.lag;
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  
  this->reciprocal_schedule_scale = another.reciprocal_schedule_scale;
  
  /*
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimators.reserve(another.measurement_covariance_estimators.size());
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i=another.measurement_covariance_estimators.begin();
       i!=another.measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimators.push_back((*i)->duplicate());
    else
      this->measurement_covariance_estimators.push_back(NULL);
  }
  */
  
  if (another.proposal!=NULL)
    this->proposal = another.proposal->independent_proposal_kernel_duplicate();
  else
    this->proposal = NULL;
  
  if (another.ensemble_factors!=NULL)
    this->ensemble_factors = another.ensemble_factors->duplicate();
  else
    this->ensemble_factors = NULL;
  
  if (another.ensemble_shifter!=NULL)
    this->ensemble_shifter = another.ensemble_shifter->duplicate();
  else
    this->ensemble_shifter = NULL;
  
  this->transform = another.transform;
  //if (another.transform!=NULL)
  //  this->transform = another.transform->duplicate();
  //else
  //  this->transform = NULL;
  
  this->packing_instructions = another.packing_instructions;
  //this->likelihood_is_evaluated = another.likelihood_is_evaluated;
  this->vector_variables = another.vector_variables;
  this->any_variables = another.any_variables;
  this->vector_variable_sizes = another.vector_variable_sizes;
  this->results_name = another.results_name;
  
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
  
  this->proposed_particles_inputted = another.proposed_particles_inputted;
  this->initial_ensemble = another.initial_ensemble;
  
  this->initialised = another.initialised;
  
  //this->transform = another.transform;
  //this->inverse_transform = another.inverse_transform;
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
  /*
  Parameters dummy_parameters = this->proposal->independent_simulate(*this->rng,
                                                                     conditioned_on_parameters);
  this->vector_variables = dummy_parameters.get_vector_variables();
  this->any_variables = dummy_parameters.get_any_variables();
  */
  
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
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
  //Ensemble* next_ensemble = current_state->add_ensemble();
  
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_ensemble,
                               index);
    //Particles particles(this->the_worker->get_particles());
    //current_state->initialise_normalised_log_weights();
    //current_state->initialise_unnormalised_log_incremental_weights(this->the_worker->get_unnormalised_log_incremental_weights());
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble,this->ensemble_factors);
  }

  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  //Ensemble* next_ensemble = current_state->add_ensemble();
  Ensemble* next_ensemble = current_state->add_ensemble(this->ensemble_factors);
  if (!this->proposed_particles_inputted)
  {
    this->the_worker->simulate(next_ensemble,
                               index,
                               conditioned_on_parameters);
  }
  else
  {
    next_ensemble->setup(this->initial_ensemble, this->ensemble_factors);
  }
  current_state->back().log_normalising_constant_ratio = 0.0;
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng),
                        this->ensemble_factors);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng)),
                        this->ensemble_factors);
  }
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  
  //EnsembleFactorVariables* simulated_ensemble_factor_variables = this->ensemble_factors->simulate_ensemble_factor_variables(simulated_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors);
}

void EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                              Particle* new_particle,
                                              const Parameters &conditioned_on_parameters) const
{
  if (this->transform==NULL)
  {
    new_particle->setup(this->proposal->independent_simulate(rng,
                                                             conditioned_on_parameters),
                        this->ensemble_factors,
                        conditioned_on_parameters);
  }
  else
  {
    new_particle->setup(this->transform->transform(this->proposal->independent_simulate(rng,
                                                                                        conditioned_on_parameters)),
                        this->ensemble_factors,
                        conditioned_on_parameters);
  }
  
  //Parameters simulated_parameters = this->proposal->independent_simulate(rng,
  //                                                                       conditioned_on_parameters);
  
  //EnsembleFactorVariables* simulated_ensemble_factor_variables = this->ensemble_factors->simulate_ensemble_factor_variables(simulated_parameters, conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  //return Particle(simulated_parameters, this->ensemble_factors, conditioned_on_parameters);
}

void EnsembleKalman::find_measurement_covariances(EnsembleKalmanOutput* simulation)
{
  simulation->back().find_measurement_covariances();
}

void EnsembleKalman::set_reciprocal_schedule_scale(double reciprocal_schedule_scale_in)
{
  this->reciprocal_schedule_scale = reciprocal_schedule_scale_in;
}

/*
void EnsembleKalman::setup_variables_using_candidate_parameters(const Parameters &candidate_parameters)
{
  this->vector_variables = candidate_parameters.get_vector_variables();
  this->packing_instructions.set_info(candidate_parameters,this->vector_variables);
  this->any_variables = candidate_parameters.get_any_variables();
}
*/
