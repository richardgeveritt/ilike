#include <numeric>
#include "ensemble_kalman_mfds.h"
#include "ensemble_kalman_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
//#include "annealed_likelihood_estimator.h"
#include "utils.h"
//#include "parameter_particle_simulator.h"
//#include "rcppparallel_smc_worker.h"
#include "sequential_ensemble_kalman_worker.h"
#include "mcmc.h"
#include "single_point_move_output.h"
#include "ensemble_factor_variables.h"
#include "generic_measurement_covariance_estimator.h"
#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "mixed_generic_direct_gaussian_measurement_covariance_estimator.h"
#include "vector_index.h"
#include "custom_distribution_proposal_kernel.h"
#include "vector_ensemble_factors.h"
#include "positive_smc_criterion.h"
//#include "standard_mcmc_output.h"
//#include "custom_distribution_proposal_kernel.h""

EnsembleKalmanMFDS::EnsembleKalmanMFDS()
   :EnsembleKalman()
{
  this->mcmc = NULL;
  this->index = NULL;
}

EnsembleKalmanMFDS::EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                                       size_t* seed_in,
                                       Data* data_in,
                                       size_t lag_in,
                                       const std::vector<Parameters> &initial_points_in,
                                       EnsembleShifter* shifter_in,
                                       double delta_t_in,
                                       size_t number_of_iterations_in,
                                       const Parameters &prior_means,
                                       const Parameters &prior_covariances,
                                       SimulateModelPtr simulate_model_in,
                                       size_t update_type,
                                       std::shared_ptr<Transform> summary_statistics_in,
                                       std::shared_ptr<Transform> transform_in,
                                       const std::vector<std::string> &measurement_variables_in,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
:EnsembleKalman(rng_in,
                seed_in,
                data_in,
                initial_points_in.size(),
                lag_in,
                shifter_in,
                transform_in,
                true,
                true,
                results_name_in), state_mean_as_data(prior_means)
{
  this->delta_t = delta_t_in;
  this->number_of_iterations = number_of_iterations_in;
  if (initial_points_in.size()==0)
    Rcpp::stop("EnsembleKalmanMFDS::EnsembleKalmanMFDS - must have at least one initial point.");
  //Parameters candidate_parameters = initial_points_in[0];
  //this->setup_variables_using_candidate_parameters(candidate_parameters);
  
  this->proposal = NULL;
  
  this->initial_ensemble = initial_points_in;
  this->proposed_particles_inputted = true;
  
  this->set_packing_instructions();
  
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  std::vector<size_t> indices;
  
  if (update_type==0)
  {
    // sequential update, with likelihood first, then prior
    measurement_covariance_estimators.reserve(2);
    measurement_covariance_estimators.push_back(new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          transform_in,
                                                                                          summary_statistics_in,
                                                                                          simulate_model_in,
                                                                                          measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    measurement_covariance_estimators.push_back(new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                          seed_in,
                                                                                                          &this->state_mean_as_data,
                                                                                                          NULL,
                                                                                                          NULL,
                                                                                                          NULL,
                                                                                                          prior_covariances.get_matrices(this->packing_instructions.states_names),
                                                                                                          this->packing_instructions.states_names));
    measurement_covariance_estimators.back()->change_data();
    
    indices.push_back(0);
    indices.push_back(1);
  }
  else if (update_type==1)
  {
    // sequential update, with prior first, then likelihood
    measurement_covariance_estimators.reserve(2);
    measurement_covariance_estimators.push_back(new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                          seed_in,
                                                                                                          &this->state_mean_as_data,
                                                                                                          NULL,
                                                                                                          NULL,
                                                                                                          NULL,
                                                                                                          prior_covariances.get_matrices(this->packing_instructions.states_names),
                                                                                                          this->packing_instructions.states_names));
    measurement_covariance_estimators.back()->change_data();
    
    measurement_covariance_estimators.push_back(new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          transform_in,
                                                                                          summary_statistics_in,
                                                                                          simulate_model_in,
                                                                                          measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    
    indices.push_back(0);
    indices.push_back(1);
  }
  else if (update_type==2)
  {
    // updating both at once, using Gaussiam for prior and intractable for llhd
    measurement_covariance_estimators.reserve(1);
    measurement_covariance_estimators.push_back(new MixedGenericDirectGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                             seed_in,
                                                                                                             data_in,
                                                                                                             transform_in,
                                                                                                             summary_statistics_in,
                                                                                                             &this->state_mean_as_data,
                                                                                                             simulate_model_in,
                                                                                                             this->packing_instructions.states_names,
                                                                                                             prior_covariances.get_matrices(this->packing_instructions.states_names),
                                                                                                             measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    indices.push_back(0);
    
  }
  else
  {
    Rcpp::stop("EnsembleKalmanMFDS::EnsembleKalmanMFDS - invalid option for sequential_update; must be 0, 1 or 2.");
  }
  
  this->index = new VectorIndex(indices);
  
  /*
  for (auto i=measurement_covariance_estimators.begin();
       i!=measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */
  
  this->ensemble_factors = new VectorEnsembleFactors(measurement_covariance_estimators);
  
  Rcpp::stop("EMFDS - need alternative to set temperature.");
  //this->ensemble_factors->set_temperature(this->delta_t);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  //this->particle_simulator = new ParameterParticleSimulator(proposal,
  //                                                          likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,
    //this->model_and_algorithm.particle_simulator,
    //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialEnsembleKalmanWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new PositiveSMCCriterion();
  this->sequencer = EnsembleSequencer(this->the_worker,
                                      schedule_in,
                                      "",
                                      25,
                                      smc_criterion);
  
  this->mcmc = NULL;
}

EnsembleKalmanMFDS::EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                                       size_t* seed_in,
                                       Data* data_in,
                                       size_t lag_in,
                                       const std::vector<Parameters> &initial_points_in,
                                       EnsembleShifter* shifter_in,
                                       double delta_t_in,
                                       size_t number_of_iterations_in,
                                       std::shared_ptr<Transform> measurement_transform_function_in,
                                       const std::vector<std::string> &measurement_variables,
                                       const std::vector<GetMatrixPtr> &measurement_noise_functions,
                                       std::shared_ptr<Transform> summary_statistics_in,
                                       std::shared_ptr<Transform> transform_in,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
:EnsembleKalman(rng_in,
                seed_in,
                data_in,
                initial_points_in.size(),
                lag_in,
                shifter_in,
                transform_in,
                true,
                true,
                results_name_in)
{
  this->delta_t = delta_t_in;
  this->number_of_iterations = number_of_iterations_in;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("EnsembleKalmanMFDS::EnsembleKalmanMFDS - must have at least one initial point.");
  //Parameters candidate_parameters = initial_points_in[0];
  //this->setup_variables_using_candidate_parameters(candidate_parameters);
  //Data candidate_measurement = transform_function_in(candidate_parameters);
  
  this->proposal = NULL;
  
  this->initial_ensemble = initial_points_in;
  this->proposed_particles_inputted = true;
  
  this->set_packing_instructions();
  
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  std::vector<size_t> indices;
  
  measurement_covariance_estimators.push_back(new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                        seed_in,
                                                                                                        data_in,
                                                                                                        transform_in,
                                                                                                        summary_statistics_in,
                                                                                                        measurement_transform_function_in,
                                                                                                        measurement_noise_functions,
                                                                                                        measurement_variables));
  measurement_covariance_estimators.back()->change_data();
    
  indices.push_back(0);

  this->index = new VectorIndex(indices);
  
  /*
  for (auto i=measurement_covariance_estimators.begin();
       i!=measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */
  
  this->ensemble_factors = new VectorEnsembleFactors(measurement_covariance_estimators);
  
  Rcpp::stop("EMFDS - need alternative to set temperature.");
  //this->ensemble_factors->set_temperature(this->delta_t);
  
  // Need to construct LikelihoodEstimator to read in to this constructor.
  //this->particle_simulator = new ParameterParticleSimulator(proposal,
  //                                                          likelihood_estimators);
  
  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,
    //this->model_and_algorithm.particle_simulator,
    //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialEnsembleKalmanWorker(this);
  }
  
  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new PositiveSMCCriterion();
  this->sequencer = EnsembleSequencer(this->the_worker,
                                      schedule_in,
                                      "",
                                      25,
                                      smc_criterion);
  
  this->mcmc = NULL;
}

EnsembleKalmanMFDS::EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                                       size_t* seed_in,
                                       Data* data_in,
                                       size_t lag_in,
                                       const std::vector<Parameters> &initial_points_in,
                                       EnsembleShifter* shifter_in,
                                       double delta_t_in,
                                       size_t number_of_iterations_in,
                                       SimulateModelPtr simulate_model_in,
                                       std::shared_ptr<Transform> summary_statistics_in,
                                       std::shared_ptr<Transform> transform_in,
                                       const std::vector<std::string> &measurement_variables_in,
                                       bool parallel_in,
                                       size_t grain_size_in,
                                       const std::string &results_name_in)
:EnsembleKalman(rng_in,
                seed_in,
                data_in,
                initial_points_in.size(),
                lag_in,
                shifter_in,
                transform_in,
                true,
                true,
                results_name_in)
{
  this->delta_t = delta_t_in;
  this->number_of_iterations = number_of_iterations_in;
  
  if (initial_points_in.size()==0)
    Rcpp::stop("EnsembleKalmanMFDS::EnsembleKalmanMFDS - must have at least one initial point.");
  //Parameters candidate_parameters = initial_points_in[0];
  //this->setup_variables_using_candidate_parameters(candidate_parameters);
  //Data candidate_measurement = transform_function_in(candidate_parameters);
  
  this->proposal = NULL;
  
  this->initial_ensemble = initial_points_in;
  this->proposed_particles_inputted = true;
  
  this->set_packing_instructions();
  
  std::vector<MeasurementCovarianceEstimator*> measurement_covariance_estimators;
  std::vector<size_t> indices;
  
  measurement_covariance_estimators.push_back(new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                        seed_in,
                                                                                        data_in,
                                                                                        transform_in,
                                                                                        summary_statistics_in,
                                                                                        simulate_model_in,
                                                                                        measurement_variables_in));
  measurement_covariance_estimators.back()->change_data();

  indices.push_back(0);

  this->index = new VectorIndex(indices);
  
  /*
  for (auto i=measurement_covariance_estimators.begin();
       i!=measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup(candidate_parameters);
  }
  */

  this->ensemble_factors = new VectorEnsembleFactors(measurement_covariance_estimators);
  
  Rcpp::stop("EMFDS - need alternative to set temperature.");
  //this->ensemble_factors->set_temperature(this->delta_t);

  if (parallel_in==TRUE)
  {
    //this->the_worker = new RcppParallelSMCWorker(this,
    //this->model_and_algorithm.particle_simulator,
    //grain_size_in);
  }
  else
  {
    this->the_worker = new SequentialEnsembleKalmanWorker(this);
  }

  std::vector<double> schedule_in;
  schedule_in.push_back(0.0);
  schedule_in.push_back(1.0);
  SMCCriterion* smc_criterion = new PositiveSMCCriterion();
  this->sequencer = EnsembleSequencer(this->the_worker,
                                      schedule_in,
                                      "",
                                      25,
                                      smc_criterion);

  this->mcmc = NULL;
}

/*
void EnsembleKalmanMFDS::setup_variables()
{
  // intentionally blank
}
*/

//Copy constructor for the EnsembleKalmanMFDS class.
EnsembleKalmanMFDS::EnsembleKalmanMFDS(const EnsembleKalmanMFDS &another)
  :EnsembleKalman(another)
{
  this->make_copy(another);
}

//Destructor for the EnsembleKalmanMFDS class.
EnsembleKalmanMFDS::~EnsembleKalmanMFDS()
{
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;
}

void EnsembleKalmanMFDS::operator=(const EnsembleKalmanMFDS &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->mcmc!=NULL)
    delete this->mcmc;
  
  if (this->index!=NULL)
    delete this->index;

  EnsembleKalman::operator=(another);
  this->make_copy(another);
}

EnsembleKalman* EnsembleKalmanMFDS::ensemble_kalman_duplicate() const
{
  return( new EnsembleKalmanMFDS(*this));
}

LikelihoodEstimator* EnsembleKalmanMFDS::duplicate() const
{
  return( new EnsembleKalmanMFDS(*this));
}

void EnsembleKalmanMFDS::make_copy(const EnsembleKalmanMFDS &another)
{
  if (another.mcmc!=NULL)
    this->mcmc = another.mcmc->mcmc_duplicate();
  else
    this->mcmc = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
  
  this->delta_t = another.delta_t;
  this->number_of_iterations = another.number_of_iterations;
  this->state_mean_as_data = another.state_mean_as_data;
  //this->sequencer = another.sequencer;
  //this->sequencer_limit_is_fixed = another.sequencer_limit_is_fixed;
}

EnsembleKalmanOutput* EnsembleKalmanMFDS::specific_run()
{
  EnsembleKalmanOutput* simulation = this->ensemble_kalman_initialise();
  this->ensemble_kalman_simulate(simulation);
  this->ensemble_kalman_evaluate(simulation);
  return simulation;
}

EnsembleKalmanOutput* EnsembleKalmanMFDS::specific_ensemble_kalman_initialise()
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this, this->lag, this->transform, this->results_name);
  return output;
}

void EnsembleKalmanMFDS::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state)
{
  if (current_state->number_of_ensemble_kalman_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            this->index);
  }
  else
  {
    // update MCMC proposals
    if (this->mcmc!=NULL)
      this->mcmc->ensemble_adapt(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    
    current_state->back().members = std::move(current_state->back().predicted_members);
    
    Ensemble* current_particles = &current_state->back();
    Ensemble* next_particles = current_state->add_ensemble(this->ensemble_factors);
    this->the_worker->move(next_particles,
                           current_particles);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
  // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
  // move in SMC also?
  //current_state->back().set_previous_target_evaluated_to_target_evaluated();
}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state)
{
  this->ensemble_kalman_evaluate_smcfixed_part(current_state);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state);
}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state)
{
  //this->the_worker->smcfixed_weight(current_state->all_particles.back());
  //current_state->initialise_next_step();
}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state)
{
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    this->the_worker->shift(&current_state->back(),
                            1.0/this->delta_t);
    
    Rcpp::stop("EMFDS - need alternative to set temperature.");
    //this->ensemble_factors->set_temperature(this->ensemble_factors->get_temperature() + this->delta_t);
    this->predict(current_state);
    this->the_worker->unpack_with_predicted(&current_state->back());
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
      current_state->write(results_name);
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    //this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    //this->sequencer.find_next_target_bisection(current_state);
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->number_of_iterations-1==current_state->enk_iteration)
    {
      //terminate = TRUE;
      current_state->terminate();
      break;
    }
    
    this->ensemble_kalman_simulate(current_state);
  }
}

MoveOutput* EnsembleKalmanMFDS::move(RandomNumberGenerator &rng,
                                     Particle &particle)
{
  
  /*
  // Outputs are created here, with memory managed by Particle hereafter.
  return EnsembleMember(simulated_parameters, simulated_ensemble_factor_variables);
  
  EnsembleMember new_member = this->mcmc->run(rng,
                                              particle);
  
  // need to delete the old factor variables if not done already
  new_member.ensemble_factor_variables = factors->simulate_ensemble_factor_variables(rng,
                                                                                     new_member.parameters);
  */
  if (this->mcmc!=NULL)
  {
    MoveOutput* mcmc_moved = this->mcmc->run(rng,
                                             particle);
    
    //mcmc_moved->back().simulate_ensemble_factor_variables(&particle);
    return mcmc_moved;
  }
  else
  {
    Particle next_particle = particle.copy_without_factor_variables();
    const EnsembleFactors* ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (ensemble_factors!=NULL)
      next_particle.simulate_ensemble_factor_variables(ensemble_factors);
    return new SinglePointMoveOutput(std::move(next_particle));
  }
  
  // also sim meas
}

/*
void EnsembleKalmanMFDS::weight_for_adapting_sequence(Ensemble &current_particles,
                                                                    double incremental_temperature)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(current_particles,
                                                      incremental_temperature);
}
*/

void EnsembleKalmanMFDS::predict(EnsembleKalmanOutput* simulation)
{
  if (this->number_of_ensemble_members>0)
  {
    simulation->back().partially_packed_predicted_members_col.clear();
    
    simulation->back().partially_packed_predicted_members_col.reserve(this->number_of_ensemble_members);
    
    arma::colvec mean = simulation->back().partially_packed_members_col[0];
    if (this->number_of_ensemble_members>1)
    {
      for (size_t i=1; i<this->number_of_ensemble_members; ++i)
      {
        mean = mean + simulation->back().partially_packed_members_col[i];
      }
    }
    mean = mean/double(this->number_of_ensemble_members);

    for (size_t i=0; i<this->number_of_ensemble_members; ++i)
    {
      //mean + sqrt(1.0/(1.0-this->delta_t))*(simulation->back().partially_packed_members_col[i]-mean);
      simulation->back().partially_packed_predicted_members_col.push_back(mean + sqrt(1.0/(1.0-this->delta_t))*(simulation->back().partially_packed_members_col[i]-mean));
    }
  }
}

EnsembleKalmanOutput* EnsembleKalmanMFDS::specific_run(const Parameters &conditioned_on_parameters)
{
  EnsembleKalmanOutput* simulation = this->ensemble_kalman_initialise(conditioned_on_parameters);
  this->ensemble_kalman_simulate(simulation,
                                 conditioned_on_parameters);
  this->ensemble_kalman_evaluate(simulation,
                                 conditioned_on_parameters);
  return simulation;
}

EnsembleKalmanOutput* EnsembleKalmanMFDS::specific_ensemble_kalman_initialise(const Parameters &conditioned_on_parameters)
{
  EnsembleKalmanOutput* output = new EnsembleKalmanOutput(this, this->lag, this->transform, this->results_name);
  return output;
}

void EnsembleKalmanMFDS::ensemble_kalman_simulate(EnsembleKalmanOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_ensemble_kalman_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            this->index,
                            conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    if (this->mcmc!=NULL)
      this->mcmc->ensemble_adapt(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_particles = &current_state->back();
    Ensemble* next_particles = current_state->add_ensemble(this->ensemble_factors);
    //this->the_worker->move(next_particles,
    //                       current_particles,
    //                       conditioned_on_parameters);
    this->the_worker->move(next_particles,
                           current_particles);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
  // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
  // move in SMC also?
  //current_state->back().set_previous_target_evaluated_to_target_evaluated();

}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate(EnsembleKalmanOutput* current_state,
                                                  const Parameters &conditioned_on_parameters)
{
  this->ensemble_kalman_evaluate_smcfixed_part(current_state,
                                               conditioned_on_parameters);
  this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                 conditioned_on_parameters);
}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                const Parameters &conditioned_on_parameters)
{
  //this->the_worker->smcfixed_weight(current_state->all_particles.back(),
  //                                  conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void EnsembleKalmanMFDS::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
                                                                                  const Parameters &conditioned_on_parameters)
{
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    this->the_worker->shift(&current_state->back(),
                            1.0/this->delta_t);
    this->predict(current_state);
    this->the_worker->unpack_with_predicted(&current_state->back());
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
      current_state->write(results_name);
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    //this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    //this->sequencer.find_next_target_bisection(current_state);
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->number_of_iterations-1==current_state->enk_iteration)
    {
      //terminate = TRUE;
      current_state->terminate();
      break;
    }
    
    this->ensemble_kalman_simulate(current_state,
                                   conditioned_on_parameters);
  }
}

void EnsembleKalmanMFDS::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* current_state)
{
  if (current_state->number_of_ensemble_kalman_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            this->index);
  }
  else
  {
    // update MCMC proposals
    if (this->mcmc!=NULL)
      this->mcmc->ensemble_adapt(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_particles = &current_state->back();
    Ensemble* next_particles = current_state->add_ensemble(this->ensemble_factors);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
  // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
  // move in SMC also?
  //current_state->back().set_previous_target_evaluated_to_target_evaluated();
  
}

void EnsembleKalmanMFDS::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* current_state,
                                                            const Parameters &conditioned_on_parameters)
{
  if (current_state->number_of_ensemble_kalman_iterations()==0)
  {
    // Simulate from the proposal.
    this->simulate_proposal(current_state,
                            this->index,
                            conditioned_on_parameters);
  }
  else
  {
    // update MCMC proposals
    if (this->mcmc!=NULL)
      this->mcmc->ensemble_adapt(current_state);
    
    // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
    Ensemble* current_particles = &current_state->back();
    Ensemble* next_particles = current_state->add_ensemble(this->ensemble_factors);
    this->the_worker->subsample_move(next_particles,
                                     current_particles);
    
    current_state->increment_enk_iteration();
    // involves complete evaluation of weights using current adaptive param
  }
  
  // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
  // move in SMC also?
  //current_state->back().set_previous_target_evaluated_to_target_evaluated();
  
}

void EnsembleKalmanMFDS::ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* current_state,
                                                            const Parameters &conditioned_on_parameters)
{
  this->ensemble_kalman_subsample_evaluate_smcfixed_part(current_state,
                                                         conditioned_on_parameters);
  this->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                           conditioned_on_parameters);
}

void EnsembleKalmanMFDS::ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* current_state,
                                                                          const Parameters &conditioned_on_parameters)
{
  //this->the_worker->subsample_smcfixed_weight(current_state->all_particles.back(),
  //                                            conditioned_on_parameters);
  //current_state->initialise_next_step();
}

void EnsembleKalmanMFDS::ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* current_state,
                                                                                            const Parameters &conditioned_on_parameters)
{
  // iterate until stop.
  bool terminate = FALSE;
  while (terminate==FALSE)
  {
    this->the_worker->pack(&current_state->back());
    this->find_measurement_covariances(current_state);
    current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_latest_log_normalising_constant_ratio();
    this->the_worker->shift(&current_state->back(),
                            1.0/this->delta_t);
    this->predict(current_state);
    this->the_worker->unpack_with_predicted(&current_state->back());
    
    current_state->set_llhd(current_state->log_likelihood);
    current_state->set_time();
    
    if (current_state->results_name!="")
      current_state->write(results_name);
    
    current_state->start_time = std::chrono::high_resolution_clock::now();
    
    //this->sequencer.find_desired_criterion(current_state);
    
    // (involves evaluating adaptive weights, using Sequencer)
    //this->sequencer.find_next_target_bisection(current_state);
    
    //this->the_worker->smcadaptive_given_smcfixed_weight(conditioned_on_parameters);
    //current_state->update_weights(this->the_worker->get_unnormalised_log_incremental_weights());
    
    // check termination, using sequencer
    if (this->number_of_iterations-1==current_state->enk_iteration)
    {
      //terminate = TRUE;
      current_state->terminate();
      break;
    }
    
    this->ensemble_kalman_simulate(current_state,
                                   conditioned_on_parameters);
  }
}

/*
MoveOutput* EnsembleKalmanMFDS::move(RandomNumberGenerator &rng,
                                     Particle &particle,
                                     const Parameters &conditioned_on_parameters)
{
  
  if (this->mcmc!=NULL)
  {
    MoveOutput* mcmc_moved = this->mcmc->run(rng,
                                             particle,
                                             conditioned_on_parameters);
    
    //mcmc_moved->back().simulate_ensemble_factor_variables(&particle,conditioned_on_parameters);
    return mcmc_moved;
  }
  else
  {
    Particle next_particle = particle.copy_without_factor_variables();
    EnsembleFactors* ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (ensemble_factors!=NULL)
      next_particle.simulate_ensemble_factor_variables(ensemble_factors,
                                                  conditioned_on_parameters);
    return new SinglePointMoveOutput(std::move(next_particle));
  }
}
*/

/*
void EnsembleKalmanMFDS::weight_for_adapting_sequence(Particles &current_particles,
                                                                    const Parameters &conditioned_on_parameters)
{
  this->the_worker->smcadaptive_given_smcfixed_weight(current_particles,
                                                      conditioned_on_parameters);
}
*/

MoveOutput* EnsembleKalmanMFDS::subsample_move(RandomNumberGenerator &rng,
                                               Particle &particle)
{
  if (this->mcmc!=NULL)
  {
    MoveOutput* mcmc_moved = this->mcmc->subsample_run(rng,
                                                       particle);
    
    //mcmc_moved->back().simulate_ensemble_factor_variables(&particle,conditioned_on_parameters);
    return mcmc_moved;
  }
  else
  {
    Particle next_particle = particle.copy_without_factor_variables();
    const EnsembleFactors* ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (ensemble_factors!=NULL)
      next_particle.subsample_simulate_ensemble_factor_variables(ensemble_factors);
    return new SinglePointMoveOutput(std::move(next_particle));
  }
}

/*
MoveOutput* EnsembleKalmanMFDS::subsample_move(RandomNumberGenerator &rng,
                                               Particle &particle,
                                               const Parameters &conditioned_on_parameters)
{
  if (this->mcmc!=NULL)
  {
    MoveOutput* mcmc_moved = this->mcmc->subsample_run(rng,
                                                       particle,
                                                       conditioned_on_parameters);
    
    //mcmc_moved->back().simulate_ensemble_factor_variables(&particle,conditioned_on_parameters);
    return mcmc_moved;
  }
  else
  {
    Particle next_particle = particle.copy_without_factor_variables();
    EnsembleFactors* ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
    if (ensemble_factors!=NULL)
      next_particle.subsample_simulate_ensemble_factor_variables(ensemble_factors,
                                                            conditioned_on_parameters);
    return new SinglePointMoveOutput(std::move(next_particle));
  }
}
*/

/*
void EnsembleKalmanMFDS::subsample_weight_for_adapting_sequence(Particles &current_particles,
                                               const Parameters &conditioned_on_parameters)
{
  this->the_worker->subsample_smcadaptive_given_smcfixed_weight(current_particles,
                                                      conditioned_on_parameters);
}
*/

// void EnsembleKalmanMFDS::smc_step(void)
// {
//
// }
//
// void EnsembleKalmanMFDS::weight_update(void)
// {
//
// }

//void EnsembleKalmanMFDS::smc_update(EnsembleKalmanOutput* current_state)
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
