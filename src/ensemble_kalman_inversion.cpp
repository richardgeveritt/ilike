#include <numeric>
#include "ensemble_kalman_inversion.h"
#include "ensemble_kalman_output.h"
#include "cess_smc_criterion.h"
#include "exact_likelihood_estimator.h"
#include "utils.h"
#include "sequential_ensemble_kalman_worker.h"
#include "mcmc.h"
#include "factor_variables.h"
#include "factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble_factors.h"
#include "move_output.h"
#include "single_point_move_output.h"
#include "generic_measurement_covariance_estimator.h"
#include "vector_ensemble_factors.h"
#include "vector_index.h"
#include "custom_distribution_proposal_kernel.h"
#include "sequential_ensemble_kalman_worker.h"
#include "ensemble_sequencer.h"
#include "positive_smc_criterion.h"
#include "ess_smc_criterion.h"
#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "direct_abc_gaussian_measurement_covariance_estimator.h"
#include "transform.h"
#include "gaussian_noise_proposal_kernel.h"

namespace ilike
{
  EnsembleKalmanInversion::EnsembleKalmanInversion()
      : EnsembleKalman()
  {
    this->mcmc.reset();
    this->index.reset();
    this->estimator_type = 0;
  }

  EnsembleKalmanInversion::EnsembleKalmanInversion(RandomNumberGenerator *rng_in,
                                                   size_t *seed_in,
                                                   Data *data_in,
                                                   size_t number_of_ensemble_members_in,
                                                   size_t lag_in,
                                                   EnsembleShifter *shifter_in,
                                                   IndependentProposalKernel *prior_in,
                                                   SimulateModelPtr simulate_model_in,
                                                   std::shared_ptr<Transform> summary_statistics_in,
                                                   std::shared_ptr<Transform> transform_in,
                                                   const std::vector<std::string> &measurement_variables_in,
                                                   double significance_level_in,
                                                   size_t estimator_type_in,
                                                   bool parallel_in,
                                                   size_t grain_size_in,
                                                   const std::string &results_name_in)
      : EnsembleKalman(rng_in,
                       seed_in,
                       data_in,
                       number_of_ensemble_members_in,
                       lag_in,
                       shifter_in,
                       transform_in,
                       true,
                       true,
                       results_name_in)
  {
    this->estimator_type = estimator_type_in;

    this->significance_level = significance_level_in;
    this->proposal.reset(prior_in);
    this->proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<MeasurementCovarianceEstimator *> measurement_covariance_estimators;
    std::vector<size_t> indices;
    measurement_covariance_estimators.reserve(1);
    measurement_covariance_estimators.push_back(new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          transform_in,
                                                                                          summary_statistics_in,
                                                                                          simulate_model_in,
                                                                                          measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->ensemble_factors.reset(new VectorEnsembleFactors(measurement_covariance_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    // this->particle_simulator = new ParameterParticleSimulator(proposal,
    //                                                          likelihood_estimators);

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
    }

    std::vector<double> schedule_in;
    schedule_in.push_back(0.0);
    schedule_in.push_back(1.0);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = EnsembleSequencer(this->the_worker.get(),
                                        schedule_in,
                                        "",
                                        25,
                                        smc_criterion);

    this->mcmc.reset();

  }

  EnsembleKalmanInversion::EnsembleKalmanInversion(RandomNumberGenerator *rng_in,
                                                   size_t *seed_in,
                                                   Data *data_in,
                                                   size_t number_of_ensemble_members_in,
                                                   size_t lag_in,
                                                   EnsembleShifter *shifter_in,
                                                   double annealing_desired_cess_in,
                                                   size_t number_of_bisections_in,
                                                   const std::string &sequence_variable_in,
                                                   const std::vector<double> &schedule_in,
                                                   IndependentProposalKernel *prior_in,
                                                   SimulateModelPtr simulate_model_in,
                                                   std::shared_ptr<Transform> summary_statistics_in,
                                                   std::shared_ptr<Transform> transform_in,
                                                   const std::vector<std::string> &measurement_variables_in,
                                                   double significance_level_in,
                                                   size_t estimator_type_in,
                                                   bool parallel_in,
                                                   size_t grain_size_in,
                                                   const std::string &results_name_in)
      : EnsembleKalman(rng_in,
                       seed_in,
                       data_in,
                       number_of_ensemble_members_in,
                       lag_in,
                       shifter_in,
                       transform_in,
                       true,
                       true,
                       results_name_in)
  {
    this->estimator_type = estimator_type_in;

    this->significance_level = significance_level_in;

    this->proposal.reset(prior_in);
    this->proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<MeasurementCovarianceEstimator *> measurement_covariance_estimators;
    std::vector<size_t> indices;
    measurement_covariance_estimators.reserve(1);
    measurement_covariance_estimators.push_back(new GenericMeasurementCovarianceEstimator(rng_in,
                                                                                          seed_in,
                                                                                          data_in,
                                                                                          transform_in,
                                                                                          summary_statistics_in,
                                                                                          simulate_model_in,
                                                                                          measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->ensemble_factors.reset(new VectorEnsembleFactors(measurement_covariance_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    // this->particle_simulator = new ParameterParticleSimulator(proposal,
    //                                                          likelihood_estimators);

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
    }

    SMCCriterion *smc_criterion = new ESSSMCCriterion(annealing_desired_cess_in);
    this->sequencer = EnsembleSequencer(this->the_worker.get(),
                                        schedule_in,
                                        sequence_variable_in,
                                        number_of_bisections_in,
                                        smc_criterion);
    this->mcmc.reset();

  }

  EnsembleKalmanInversion::EnsembleKalmanInversion(RandomNumberGenerator *rng_in,
                                                   size_t *seed_in,
                                                   Data *data_in,
                                                   size_t number_of_ensemble_members_in,
                                                   size_t lag_in,
                                                   EnsembleShifter *shifter_in,
                                                   double annealing_desired_cess_in,
                                                   size_t number_of_bisections_in,
                                                   const std::string &sequence_variable_in,
                                                   const std::vector<double> &schedule_in,
                                                   IndependentProposalKernel *prior_in,
                                                   std::shared_ptr<Transform> measurement_transform_function_in,
                                                   const std::vector<std::string> &measurement_variables,
                                                   const std::vector<GetMatrixPtr> &measurement_noise_functions_in,
                                                   std::shared_ptr<Transform> summary_statistics_in,
                                                   std::shared_ptr<Transform> transform_in,
                                                   double significance_level_in,
                                                   size_t estimator_type_in,
                                                   bool parallel_in,
                                                   size_t grain_size_in,
                                                   const std::string &results_name_in)
      : EnsembleKalman(rng_in,
                       seed_in,
                       data_in,
                       number_of_ensemble_members_in,
                       lag_in,
                       shifter_in,
                       transform_in,
                       true,
                       true,
                       results_name_in)
  {
    this->estimator_type = estimator_type_in;

    this->significance_level = significance_level_in;

    this->proposal.reset(prior_in);
    this->proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<MeasurementCovarianceEstimator *> measurement_covariance_estimators;
    std::vector<size_t> indices;
    measurement_covariance_estimators.reserve(1);
    measurement_covariance_estimators.push_back(new DirectNonLinearGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                          seed_in,
                                                                                                          data_in,
                                                                                                          transform_in,
                                                                                                          summary_statistics_in,
                                                                                                          measurement_transform_function_in,
                                                                                                          measurement_noise_functions_in,
                                                                                                          measurement_variables));
    measurement_covariance_estimators.back()->change_data();
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->ensemble_factors.reset(new VectorEnsembleFactors(measurement_covariance_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    // this->particle_simulator = new ParameterParticleSimulator(proposal,
    //                                                          likelihood_estimators);

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
    }

    SMCCriterion *smc_criterion = new ESSSMCCriterion(annealing_desired_cess_in);
    this->sequencer = EnsembleSequencer(this->the_worker.get(),
                                        schedule_in,
                                        sequence_variable_in,
                                        number_of_bisections_in,
                                        smc_criterion);
    this->mcmc.reset();

  }

  EnsembleKalmanInversion::EnsembleKalmanInversion(RandomNumberGenerator *rng_in,
                                                   size_t *seed_in,
                                                   Data *data_in,
                                                   size_t number_of_ensemble_members_in,
                                                   size_t lag_in,
                                                   EnsembleShifter *shifter_in,
                                                   size_t number_of_targets,
                                                   // size_t number_of_bisections_in,
                                                   const std::string &sequence_variable_in,
                                                   IndependentProposalKernel *prior_in,
                                                   const std::vector<std::string> &measurement_variables_in,
                                                   double min_epsilon_in,
                                                   const std::string &scale_variable_in,
                                                   const arma::colvec &scale_in,
                                                   std::shared_ptr<Transform> summary_statistics_in,
                                                   std::shared_ptr<Transform> transform_in,
                                                   double significance_level_in,
                                                   size_t estimator_type_in,
                                                   bool parallel_in,
                                                   size_t grain_size_in,
                                                   const std::string &results_name_in)
      : EnsembleKalman(rng_in,
                       seed_in,
                       data_in,
                       number_of_ensemble_members_in,
                       lag_in,
                       shifter_in,
                       transform_in,
                       true,
                       true,
                       results_name_in)
  {
    this->estimator_type = estimator_type_in;

    this->significance_level = significance_level_in;

    this->proposal.reset(prior_in);
    this->proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<MeasurementCovarianceEstimator *> measurement_covariance_estimators;
    std::vector<size_t> indices;
    measurement_covariance_estimators.reserve(1);
    measurement_covariance_estimators.push_back(new DirectABCGaussianMeasurementCovarianceEstimator(rng_in,
                                                                                                    seed_in,
                                                                                                    data_in,
                                                                                                    transform_in,
                                                                                                    summary_statistics_in,
                                                                                                    min_epsilon_in,
                                                                                                    // sequence_variable_in,
                                                                                                    scale_variable_in,
                                                                                                    measurement_variables_in));
    measurement_covariance_estimators.back()->change_data();
    indices.push_back(0);
    this->index.reset(new VectorIndex(indices));

    this->ensemble_factors.reset(new VectorEnsembleFactors(measurement_covariance_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    // this->particle_simulator = new ParameterParticleSimulator(proposal,
    //                                                          likelihood_estimators);

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
    }

    // using the vector to store the information needed later, when we can set up the sequence
    std::vector<double> schedule_in;
    schedule_in.reserve(number_of_targets);
    for (size_t i = 0; i < number_of_targets + 1; ++i)
      schedule_in.push_back(min_epsilon_in);
    SMCCriterion *smc_criterion = new PositiveSMCCriterion();
    this->sequencer = EnsembleSequencer(this->the_worker.get(),
                                        schedule_in,
                                        sequence_variable_in,
                                        100, // number_of_bisections_in,
                                        smc_criterion);

    // sort out ABC scaling

    size_t data_size = 0;
    for (auto i = measurement_variables_in.begin();
         i != measurement_variables_in.end();
         ++i)
    {
      data_size = data_size + (*this->data)[*i].n_elem;
    }

    // if scale variable is not set, we are not using a scale at all, and we set it to be ones
    if (scale_variable_in == "")
    {
      this->sequencer.scale_variable = "default_scale";
      arma::colvec scale = arma::colvec(data_size);
      scale.fill(1.0);
      this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale;
      this->sequencer.find_scale = false;
    }
    else // if the scale variable is set, then either we have read in a valid scale, or it must be estimated through simulation
    {
      this->sequencer.scale_variable = scale_variable_in;

      if (scale_in.n_elem == data_size)
      {
        this->sequencer.schedule_parameters[this->sequencer.scale_variable] = scale_in;
        this->sequencer.find_scale = false;
      }
      else
      {
        Rcpp::stop("The scale variable must be the same length as the data.");
      }
    }

    this->mcmc.reset();

    this->set_reciprocal_schedule_scale(min_epsilon_in);

  }

  EnsembleKalmanInversion::EnsembleKalmanInversion(RandomNumberGenerator *rng_in,
                                                   size_t *seed_in,
                                                   Data *data_in,
                                                   size_t number_of_ensemble_members_in,
                                                   size_t lag_in,
                                                   EnsembleShifter *shifter_in,
                                                   SMCCriterion *adaptive_target_criterion_in,
                                                   size_t number_of_bisections_in,
                                                   SMCTermination *termination_in,
                                                   const std::string &sequence_variable_in,
                                                   const std::vector<double> &schedule_in,
                                                   IndependentProposalKernel *prior_in,
                                                   const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                                                   const std::vector<MeasurementCovarianceEstimator *> &estimators_in,
                                                   std::shared_ptr<Transform> transform_in,
                                                   double significance_level_in,
                                                   size_t estimator_type_in,
                                                   bool parallel_in,
                                                   size_t grain_size_in,
                                                   const std::string &results_name_in)
      : EnsembleKalman(rng_in,
                       seed_in,
                       data_in,
                       number_of_ensemble_members_in,
                       lag_in,
                       shifter_in,
                       transform_in,
                       true,
                       true,
                       results_name_in)
  {
    this->estimator_type = estimator_type_in;

    this->significance_level = significance_level_in;

    this->proposal.reset(prior_in);
    this->proposal->set_proposal_parameters(&this->algorithm_parameters);

    std::vector<MeasurementCovarianceEstimator *> measurement_covariance_estimators = estimators_in;

    std::vector<size_t> indices;

    for (size_t i = 0;
         i < measurement_covariance_estimators.size();
         ++i)
    {
      measurement_covariance_estimators[i]->change_data();
      indices.push_back(i);
    }

    this->index.reset(new VectorIndex(indices));

    this->ensemble_factors.reset(new VectorEnsembleFactors(measurement_covariance_estimators));

    // Need to construct LikelihoodEstimator to read in to this constructor.
    // this->particle_simulator = new ParameterParticleSimulator(proposal,
    //                                                          likelihood_estimators);

    if (parallel_in == TRUE)
    {
    }
    else
    {
      this->the_worker.reset(new SequentialEnsembleKalmanWorker(this));
    }

    this->sequencer = EnsembleSequencer(this->the_worker.get(),
                                        schedule_in,
                                        sequence_variable_in,
                                        number_of_bisections_in,
                                        adaptive_target_criterion_in,
                                        termination_in);
    this->mcmc.reset();

  }

  // Copy constructor for the EnsembleKalmanInversion class.
  EnsembleKalmanInversion::EnsembleKalmanInversion(const EnsembleKalmanInversion &another)
      : EnsembleKalman(another)
  {
    this->make_copy(another);
  }

  // Destructor for the EnsembleKalmanInversion class.
  EnsembleKalmanInversion::~EnsembleKalmanInversion()
  {
    this->mcmc.reset();

    this->index.reset();
  }

  void EnsembleKalmanInversion::operator=(const EnsembleKalmanInversion &another)
  {
    if (this == &another)
    { // if a==a
      return;
    }

    this->mcmc.reset();

    this->index.reset();

    EnsembleKalman::operator=(another);
    this->make_copy(another);
  }

  EnsembleKalman *EnsembleKalmanInversion::ensemble_kalman_duplicate() const
  {
    return (new EnsembleKalmanInversion(*this));
  }

  LikelihoodEstimator *EnsembleKalmanInversion::duplicate() const
  {
    return (new EnsembleKalmanInversion(*this));
  }

  void EnsembleKalmanInversion::make_copy(const EnsembleKalmanInversion &another)
  {
    if (another.mcmc != NULL)
      this->mcmc.reset(another.mcmc->mcmc_duplicate());
    else
      this->mcmc.reset();

    if (another.index != NULL)
      this->index.reset(another.index->duplicate());
    else
      this->index.reset();

    this->significance_level = another.significance_level;

    this->estimator_type = another.estimator_type;
  }

  EnsembleKalmanOutput *EnsembleKalmanInversion::specific_run()
  {
    EnsembleKalmanOutput *simulation = this->ensemble_kalman_initialise();
    this->ensemble_kalman_simulate(simulation);
    this->ensemble_kalman_evaluate(simulation);
    return simulation;
  }

  EnsembleKalmanOutput *EnsembleKalmanInversion::specific_ensemble_kalman_initialise()
  {
    EnsembleKalmanOutput *output = new EnsembleKalmanOutput(this, this->lag, this->transform, this->results_name);
    return output;
  }

  void EnsembleKalmanInversion::ensemble_kalman_simulate(EnsembleKalmanOutput *current_state)
  {
    if (current_state->all_ensembles.size() == 0)
    {
      // Simulate from the proposal.
      this->simulate_proposal(current_state,
                              this->index.get());
    }
    else
    {
      // update MCMC proposals
      if (this->mcmc != NULL)
        this->mcmc->ensemble_adapt(current_state);

      // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
      Ensemble *current_particles = &current_state->back();
      Ensemble *next_particles = current_state->add_ensemble(this->ensemble_factors.get());
      this->the_worker->move(next_particles,
                             current_particles);

      current_state->increment_enk_iteration();

      // involves complete evaluation of weights using current adaptive param
    }

    // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
    // move in SMC also?
    // current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate(EnsembleKalmanOutput *current_state)
  {
    this->ensemble_kalman_evaluate_smcfixed_part(current_state);
    this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state);
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput *current_state)
  {
  }

  void EnsembleKalmanInversion::set_abc_schedule(EnsembleKalmanOutput *current_state)
  {
    std::vector<double> dummy_schedule = this->sequencer.schedule;
    size_t T = dummy_schedule.size() - 1;
    std::vector<double> abc_schedule;

    // Argument could be read in or inferred from simulated points - the same across MCMC iterations.

    arma::colvec abc_scale;
    arma::colvec simulation_scale;
    if (this->sequencer.find_scale == true)
    {
      abc_scale = arma::stddev(current_state->back().packed_members, 0, 0);
      this->sequencer.schedule_parameters[this->sequencer.scale_variable] = abc_scale;
      simulation_scale = abc_scale;
    }
    else
    {
      abc_scale = this->sequencer.schedule_parameters[this->sequencer.scale_variable];

      // Scale found from simulations from the prior - different depending on the parameter we look at.
      simulation_scale = arma::stddev(current_state->back().packed_members, 0, 0).as_col();
    }

    double epsilon = dummy_schedule.back();
    double epsilon2 = pow(epsilon, 2.0);

    for (size_t i = 0; i < abc_scale.n_elem; ++i)
    {
      if (abc_scale[i] == 0.0)
      {
        abc_scale[i] = simulation_scale[i];
      }
    }
    arma::colvec ratio = simulation_scale / abc_scale;

    double kappa = arma::mean(ratio);
    double kappa2 = pow(kappa, 2.0);

    double min_temperature = 0.0;

    if (kappa < epsilon)
    {
      abc_schedule.clear();
      abc_schedule.push_back(min_temperature);
      abc_schedule.push_back(1.0);
    }
    else
    {
      abc_schedule.reserve(T + 1);
      abc_schedule.push_back(min_temperature);
      for (size_t i = 1; i < T; ++i)
      {
        double t = double(i) / double(T);
        double alpha = exp(2.0 * log(kappa / epsilon) * t + log(epsilon2) - log(kappa2 - epsilon2)) - epsilon2 / (kappa2 - epsilon2);
        abc_schedule.push_back(alpha);
      }
      abc_schedule.push_back(1.0);
    }

    this->sequencer.set_schedule(abc_schedule);
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput *current_state)
  {
    // set sequencer to have values from conditioned_on_parameters
    if (!this->sequencer_limit_is_fixed)
      Rcpp::stop("EnsembleKalmanInversion::evaluate_smcadaptive_part_given_smcfixed_smc - need fixed sequencer limit if we are not conditioning on parameters.");

    // iterate until stop.
    bool terminate = FALSE;
    while (terminate == FALSE)
    {
      this->the_worker->pack(&current_state->back());

      // A hack so that when using in ABC, the sequence is determined using the approach described in the paper.
      // Relies on the ABC scale being found from the sample.
      double path_sampling_U_multiplier = 1.0;
      double temperature;

      if ((this->reciprocal_schedule_scale != 0.0) && (current_state->number_of_ensemble_kalman_iterations() == 1))
      {

        this->set_abc_schedule(current_state);
        temperature = this->sequencer.schedule[0]; // this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      }

      if (this->significance_level < 1.0)
      {
        current_state->skip_to_end_of_sequence_if_points_are_gaussian(this->significance_level);
      }

      this->find_measurement_covariances(current_state);
      this->sequencer.find_desired_criterion(current_state);
      // (involves evaluating adaptive weights, using Sequencer)
      this->sequencer.find_next_target_bisection(current_state,
                                                 this->index.get());

      // need to make sure that everything is called for first ensemble (for path sampling)
      // if (current_state->number_of_ensemble_kalman_iterations()<=1)
      //  int thing = 1;

      // Kalman gains are calculated as a part of this step, since (for options 1 and 2) the unconditional measurement covariance needs to be be computed here. The ESS is also computed here.
      std::vector<double> initial_log_measurement_likelihood_means;
      std::vector<double> initial_log_measurement_likelihood_variances;
      if (this->estimator_type == 0)
      {
        // the Monte Carlo estimate of the increment does not give the Kalman
        // gains that the shift below needs (without them the ensemble never
        // moves); computed first, since calculate_kalman_gains also sets the
        // ESS, which the Monte Carlo estimate then sets from its weights
        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_mc_inversion_latest_log_normalising_constant_ratio(this->index.get(),
                                                                                                                                                    1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 1)
      {
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 2)
      {
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_unbiased_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path1_inversion_initial_quantities(initial_log_measurement_likelihood_means,
                                                                      temperature,
                                                                      path_sampling_U_multiplier);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path2_inversion_initial_quantities(initial_log_measurement_likelihood_means, initial_log_measurement_likelihood_variances);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }

      this->the_worker->shift(&current_state->back(),
                              1.0 / this->sequencer.schedule_difference);
      this->the_worker->unpack(&current_state->back());

      temperature = this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      double path_sampling_inverse_incremental_temperature = 1.0 / this->sequencer.schedule_difference;

      if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                         path_sampling_inverse_incremental_temperature,
                                                                                                                                                         temperature,
                                                                                                                                                         path_sampling_U_multiplier);
        }
        else
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                         path_sampling_inverse_incremental_temperature,
                                                                                                                                                         temperature,
                                                                                                                                                         path_sampling_U_multiplier);
        }
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                         initial_log_measurement_likelihood_variances,
                                                                                                                                                         path_sampling_inverse_incremental_temperature);
        }
        else
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                         current_state->back().log_measurement_likelihood_variances,
                                                                                                                                                         path_sampling_inverse_incremental_temperature);
        }
      }

      current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

      current_state->set_llhd(current_state->log_likelihood);
      current_state->set_time();

      if (current_state->results_name != "")
      {
        current_state->write(results_name);
      }

      current_state->start_time = std::chrono::high_resolution_clock::now();

      // check termination, using sequencer
      if (this->sequencer.check_termination())
      {
        current_state->terminate();
        break;
      }

      this->ensemble_kalman_simulate(current_state);
    }
  }

  MoveOutput *EnsembleKalmanInversion::move(RandomNumberGenerator &rng,
                                            Particle &particle)
  {

    if (this->mcmc != NULL)
    {
      MoveOutput *mcmc_moved = this->mcmc->run(rng,
                                               particle);

      return mcmc_moved;
    }
    else
    {
      Particle next_particle = particle.copy_without_factor_variables();
      const EnsembleFactors *ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
      if (ensemble_factors != NULL)
        next_particle.simulate_ensemble_factor_variables(ensemble_factors);
      return new SinglePointMoveOutput(std::move(next_particle));
    }

    // also sim meas
  }

  EnsembleKalmanOutput *EnsembleKalmanInversion::specific_run(const Parameters &conditioned_on_parameters)
  {
    EnsembleKalmanOutput *simulation = this->ensemble_kalman_initialise(conditioned_on_parameters);
    this->ensemble_kalman_simulate(simulation,
                                   conditioned_on_parameters);
    this->ensemble_kalman_evaluate(simulation,
                                   conditioned_on_parameters);
    return simulation;
  }

  EnsembleKalmanOutput *EnsembleKalmanInversion::specific_ensemble_kalman_initialise(const Parameters &conditioned_on_parameters)
  {
    EnsembleKalmanOutput *output = new EnsembleKalmanOutput(this, this->lag, this->transform, this->results_name);
    return output;
  }

  void EnsembleKalmanInversion::ensemble_kalman_simulate(EnsembleKalmanOutput *current_state,
                                                         const Parameters &conditioned_on_parameters)
  {
    if (current_state->all_ensembles.size() == 0)
    {
      // Simulate from the proposal.
      this->simulate_proposal(current_state,
                              this->index.get(),
                              conditioned_on_parameters);
    }
    else
    {
      // update MCMC proposals
      if (this->mcmc != NULL)
        this->mcmc->ensemble_adapt(current_state);

      // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
      Ensemble *current_particles = &current_state->back();
      Ensemble *next_particles = current_state->add_ensemble(this->ensemble_factors.get());
      this->the_worker->move(next_particles,
                             current_particles);

      current_state->increment_enk_iteration();
      // involves complete evaluation of weights using current adaptive param
    }

    // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
    // move in SMC also?
    // current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate(EnsembleKalmanOutput *current_state,
                                                         const Parameters &conditioned_on_parameters)
  {
    this->ensemble_kalman_evaluate_smcfixed_part(current_state,
                                                 conditioned_on_parameters);
    this->ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                   conditioned_on_parameters);
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput *current_state,
                                                                       const Parameters &conditioned_on_parameters)
  {
  }

  void EnsembleKalmanInversion::ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput *current_state,
                                                                                         const Parameters &conditioned_on_parameters)
  {
    // set sequencer to have values from conditioned_on_parameters
    if (!this->sequencer_limit_is_fixed)
      this->sequencer.set_next_with_parameter(conditioned_on_parameters);
    else
      this->sequencer.reset();

    // iterate until stop.
    bool terminate = FALSE;
    while (terminate == FALSE)
    {
      this->the_worker->pack(&current_state->back());

      // A hack so that when using in ABC, the sequence is determined using the approach described in the paper.
      // Relies on the ABC scale being found from the sample.

      double path_sampling_U_multiplier = 1.0;
      double temperature;

      if ((this->reciprocal_schedule_scale != 0.0) && (current_state->number_of_ensemble_kalman_iterations() == 1))
      {
        this->set_abc_schedule(current_state);
        temperature = this->sequencer.schedule[0]; // this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      }

      if (this->significance_level < 1.0)
      {
        current_state->skip_to_end_of_sequence_if_points_are_gaussian(this->significance_level);
      }

      this->find_measurement_covariances(current_state);

      this->sequencer.find_desired_criterion(current_state);

      // (involves evaluating adaptive weights, using Sequencer)
      this->sequencer.find_next_target_bisection(current_state,
                                                 this->index.get());

      // Kalman gains are calculated as a part of this step, since (for options 1 and 2) the unconditional measurement covariance needs to be be computed here. The ESS is also computed here.
      std::vector<double> initial_log_measurement_likelihood_means;
      std::vector<double> initial_log_measurement_likelihood_variances;
      if (this->estimator_type == 0)
      {
        // (gains first: see the unconditioned version above)
        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_mc_inversion_latest_log_normalising_constant_ratio(this->index.get(),1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 1)
      {
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 2)
      {
        current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_unbiased_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path1_inversion_initial_quantities(initial_log_measurement_likelihood_means,
                                                                      temperature,
                                                                      path_sampling_U_multiplier);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path2_inversion_initial_quantities(initial_log_measurement_likelihood_means, initial_log_measurement_likelihood_variances);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }

      this->the_worker->shift(&current_state->back(),
                              1.0 / this->sequencer.schedule_difference);
      this->the_worker->unpack(&current_state->back());

      temperature = this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      double path_sampling_inverse_incremental_temperature = 1.0 / this->sequencer.schedule_difference;

      if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                         path_sampling_inverse_incremental_temperature,
                                                                                                                                                         temperature,
                                                                                                                                                         path_sampling_U_multiplier);
        }
        else
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                         path_sampling_inverse_incremental_temperature,
                                                                                                                                                         temperature,
                                                                                                                                                         path_sampling_U_multiplier);
        }
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                         initial_log_measurement_likelihood_variances,
                                                                                                                                                         path_sampling_inverse_incremental_temperature);
        }
        else
        {
          current_state->log_likelihood = current_state->log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                         current_state->back().log_measurement_likelihood_variances,
                                                                                                                                                         path_sampling_inverse_incremental_temperature);
        }
      }

      current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

      current_state->set_llhd(current_state->log_likelihood);
      current_state->set_time();

      if (current_state->results_name != "")
      {
        current_state->write(results_name);
      }
      current_state->start_time = std::chrono::high_resolution_clock::now();

      // check termination, using sequencer
      if (this->sequencer.check_termination())
      {
        current_state->terminate();
        break;
      }

      this->ensemble_kalman_simulate(current_state,
                                     conditioned_on_parameters);
    }
  }

  void EnsembleKalmanInversion::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput *current_state)
  {
    if (current_state->all_ensembles.size() == 0)
    {
      // Simulate from the proposal.
      this->simulate_proposal(current_state,
                              this->index.get());
    }
    else
    {
      // update MCMC proposals
      if (this->mcmc != NULL)
        this->mcmc->ensemble_adapt(current_state);

      // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
      Ensemble *current_particles = &current_state->back();
      Ensemble *next_particles = current_state->add_ensemble(this->ensemble_factors.get());
      this->the_worker->move(next_particles,
                             current_particles);

      current_state->increment_enk_iteration();
      // involves complete evaluation of weights using current adaptive param
    }

    // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
    // move in SMC also?
    // current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }

  void EnsembleKalmanInversion::ensemble_kalman_subsample_simulate(EnsembleKalmanOutput *current_state,
                                                                   const Parameters &conditioned_on_parameters)
  {
    if (current_state->all_ensembles.size() == 0)
    {
      // Simulate from the proposal.
      this->simulate_proposal(current_state,
                              this->index.get(),
                              conditioned_on_parameters);
    }
    else
    {
      // update MCMC proposals
      if (this->mcmc != NULL)
        this->mcmc->ensemble_adapt(current_state);

      // move (sometimes only do this when resample - to do this, adapt number of moves based on diversity of positions);
      Ensemble *current_particles = &current_state->back();
      Ensemble *next_particles = current_state->add_ensemble(this->ensemble_factors.get());
      this->the_worker->move(next_particles,
                             current_particles);

      current_state->increment_enk_iteration();
      // involves complete evaluation of weights using current adaptive param
    }

    // moved here - different to SMC due to setting target evaluated rather than previous target evaluated in first step
    // move in SMC also?
    // current_state->back().set_previous_target_evaluated_to_target_evaluated();
  }

  void EnsembleKalmanInversion::ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput *current_state,
                                                                   const Parameters &conditioned_on_parameters)
  {
    this->ensemble_kalman_subsample_evaluate_smcfixed_part(current_state,
                                                           conditioned_on_parameters);
    this->ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(current_state,
                                                                             conditioned_on_parameters);
  }

  void EnsembleKalmanInversion::ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput *current_state,
                                                                                 const Parameters &conditioned_on_parameters)
  {
  }

  void EnsembleKalmanInversion::ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput *current_state,
                                                                                                   const Parameters &conditioned_on_parameters)
  {
    if (!this->sequencer_limit_is_fixed)
      this->sequencer.set_next_with_parameter(conditioned_on_parameters);
    else
      this->sequencer.reset();

    // iterate until stop.
    bool terminate = FALSE;
    while (terminate == FALSE)
    {
      this->the_worker->pack(&current_state->back());

      // A hack so that when using in ABC, the sequence is determined using the approach described in the paper.
      // Relies on the ABC scale being found from the sample.
      double path_sampling_U_multiplier = 1.0;
      double temperature;

      if ((this->reciprocal_schedule_scale != 0.0) && (current_state->number_of_ensemble_kalman_iterations() == 1))
      {
        this->set_abc_schedule(current_state);
        temperature = this->sequencer.schedule[0]; // this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      }

      if (this->significance_level < 1.0)
      {
        current_state->skip_to_end_of_sequence_if_points_are_gaussian(this->significance_level);
      }

      this->find_measurement_covariances(current_state);

      this->sequencer.find_desired_criterion(current_state);

      // (involves evaluating adaptive weights, using Sequencer)
      this->sequencer.subsample_find_next_target_bisection(current_state,
                                                           this->index.get());

      // Kalman gains are calculated as a part of this step, since (for options 1 and 2) the unconditional measurement covariance needs to be be computed here. The ESS is also computed here.
      // Kalman gains are calculated as a part of this step, since (for options 1 and 2) the unconditional measurement covariance needs to be be computed here. The ESS is also computed here.
      std::vector<double> initial_log_measurement_likelihood_means;
      std::vector<double> initial_log_measurement_likelihood_variances;
      if (this->estimator_type == 0)
      {
        // (gains first: see the unconditioned version above)
        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
        current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_mc_inversion_latest_log_normalising_constant_ratio(this->index.get(),
                                                                                                                                                                        1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 1)
      {
        current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 2)
      {
        current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_unbiased_inversion_latest_log_normalising_constant_ratio(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path1_inversion_initial_quantities(initial_log_measurement_likelihood_means,
                                                                      temperature,
                                                                      path_sampling_U_multiplier);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->calculate_path2_inversion_initial_quantities(initial_log_measurement_likelihood_means, initial_log_measurement_likelihood_variances);
        }

        current_state->calculate_kalman_gains(1.0 / this->sequencer.schedule_difference);
      }

      this->the_worker->shift(&current_state->back(),
                              1.0 / this->sequencer.schedule_difference);
      this->the_worker->unpack(&current_state->back());

      temperature = this->sequencer.schedule_parameters[this->sequencer.variable_name][0];
      double path_sampling_inverse_incremental_temperature = 1.0 / this->sequencer.schedule_difference;

      if (this->estimator_type == 3)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                                             path_sampling_inverse_incremental_temperature,
                                                                                                                                                                             temperature,
                                                                                                                                                                             path_sampling_U_multiplier);
        }
        else
        {
          current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_path1_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                                             path_sampling_inverse_incremental_temperature,
                                                                                                                                                                             temperature,
                                                                                                                                                                             path_sampling_U_multiplier);
        }
      }
      else if (this->estimator_type == 4)
      {
        if (current_state->number_of_ensemble_kalman_iterations() == 1)
        {
          current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(initial_log_measurement_likelihood_means,
                                                                                                                                                                             initial_log_measurement_likelihood_variances,
                                                                                                                                                                             path_sampling_inverse_incremental_temperature);
        }
        else
        {
          current_state->subsample_log_likelihood = current_state->subsample_log_likelihood + current_state->calculate_path2_inversion_latest_log_normalising_constant_ratio(current_state->back().log_measurement_likelihood_means,
                                                                                                                                                                             current_state->back().log_measurement_likelihood_variances,
                                                                                                                                                                             path_sampling_inverse_incremental_temperature);
        }
      }

      current_state->back().schedule_parameters = this->sequencer.schedule_parameters.deep_copy();

      current_state->set_llhd(current_state->subsample_log_likelihood);
      current_state->set_time();

      if (current_state->results_name != "")
      {
        current_state->write(results_name);
      }
      current_state->start_time = std::chrono::high_resolution_clock::now();

      // check termination, using sequencer
      if (this->sequencer.check_termination())
      {
        current_state->terminate();
        break;
      }

      this->ensemble_kalman_subsample_simulate(current_state,
                                               conditioned_on_parameters);
    }
  }

  MoveOutput *EnsembleKalmanInversion::subsample_move(RandomNumberGenerator &rng,
                                                      Particle &particle)
  {
    if (this->mcmc != NULL)
    {
      MoveOutput *mcmc_moved = this->mcmc->subsample_run(rng,
                                                         particle);

      return mcmc_moved;
    }
    else
    {
      Particle next_particle = particle.copy_without_factor_variables();
      const EnsembleFactors *ensemble_factors = particle.ensemble_factor_variables->get_ensemble_factors();
      if (ensemble_factors != NULL)
        next_particle.subsample_simulate_ensemble_factor_variables(ensemble_factors);
      return new SinglePointMoveOutput(std::move(next_particle));
    }
  }

}
