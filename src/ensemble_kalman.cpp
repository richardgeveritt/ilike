#include "ensemble_kalman.h"
#include "ensemble_kalman_output.h"
#include "ensemble_kalman_worker.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "independent_proposal_kernel.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "ensemble_shifter.h"

EnsembleKalman::EnsembleKalman()
  :LikelihoodEstimator()
{
  this->proposal = NULL;
  this->ensemble_factors = NULL;
  this->ensemble_shifter = NULL;
  //this->measurement_covariance_estimators.resize(0);
  //this->measurements_names.resize(0);
  //this->states_names.resize(0);
}

EnsembleKalman::EnsembleKalman(RandomNumberGenerator* rng_in,
                               size_t* seed_in,
                               Data* data_in,
                               bool smcfixed_flag_in,
                               bool sequencer_limit_is_fixed_in)
:LikelihoodEstimator(rng_in, seed_in, data_in)
{
  this->proposal = NULL;
  this->ensemble_factors = NULL;
  this->ensemble_shifter = NULL;
  this->sequencer_limit_is_fixed = sequencer_limit_is_fixed_in;
  this->smcfixed_flag = smcfixed_flag_in;
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
  this->smcfixed_flag = another.smcfixed_flag;
  this->number_of_ensemble_members = another.number_of_ensemble_members;
  this->lag = another.lag;
  this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
  
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
  
  this->packing_instructions = another.packing_instructions;
  this->likelihood_is_evaluated = another.likelihood_is_evaluated;
}

EnsembleKalmanOutput* EnsembleKalman::run()
{
  return this->specific_run();
}

EnsembleKalmanOutput* EnsembleKalman::run(const Parameters &conditioned_on_parameters)
{
  return this->specific_run(conditioned_on_parameters);
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise()
{
  return this->ensemble_kalman_initialise();
}

LikelihoodEstimatorOutput* EnsembleKalman::initialise(const Parameters &conditioned_on_parameters)
{
  return this->ensemble_kalman_initialise(conditioned_on_parameters);
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index)
{
  // Simulate from proposal.
  Ensemble* next_ensemble = current_state->add_ensemble();
  this->the_worker->simulate(next_ensemble,
                             index);
}

void EnsembleKalman::simulate_proposal(EnsembleKalmanOutput* current_state,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters)
{
  // Simulate from proposal.
  Ensemble* current_ensemble = current_state->add_ensemble();
  this->the_worker->simulate(current_ensemble,
                             index,
                             conditioned_on_parameters);
}

Particle EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng) const
{
  Parameters simulated_parameters = this->proposal->independent_simulate(rng);
  
  EnsembleFactorVariables* simulated_ensemble_factor_variables = this->ensemble_factors->simulate_ensemble_factor_variables(simulated_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, simulated_ensemble_factor_variables);
}

Particle EnsembleKalman::simulate_ensemble_member(RandomNumberGenerator &rng,
                                                  const Parameters &conditioned_on_parameters) const
{
  Parameters simulated_parameters = this->proposal->independent_simulate(rng,
                                                                         conditioned_on_parameters);
  
  EnsembleFactorVariables* simulated_ensemble_factor_variables = this->ensemble_factors->simulate_ensemble_factor_variables(simulated_parameters,
                                                                                                                            conditioned_on_parameters);
  
  // Outputs are created here, with memory managed by Particle hereafter.
  return Particle(simulated_parameters, simulated_ensemble_factor_variables);
}

void EnsembleKalman::find_measurement_covariances(EnsembleKalmanOutput* simulation)
{
  simulation->back().find_measurement_covariances();
}
