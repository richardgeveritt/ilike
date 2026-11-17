#include "recipe_internal.h"
#include "recipe_from_r.h"

// R binding (dev/plan.md 3.14, I1): R lists of values for the state
// (initial values, algorithm parameters) become Parameters here; Parameters
// keeps every kind of value it supports, including "any" values.
std::vector<Parameters> make_initial_points(const Rcpp::List &initial_values) {
  std::vector<Parameters> output;
  output.reserve(initial_values.size());

  for (size_t i = 0; i < initial_values.size(); ++i) {
    Parameters current_parameters;
    if (Rf_isNewList(initial_values[i])) {
      Rcpp::List list_params = initial_values[i];
      if (list_params.size() > 0) {
        CharacterVector names = list_params.names();

        for (size_t i = 0; i < names.size(); ++i) {
          // (a matrix first: the vector test is also true for a matrix,
          // which then lost its shape)
          if (Rf_isMatrix(list_params[Rcpp::as<std::string>(names[i])])) {
            arma::mat param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          } else if (isDoubleInList(recipe_from_R(list_params), Rcpp::as<std::string>(names[i]))) {
            double param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          } else if (isNumericVectorInList(recipe_from_R(list_params),
                                           Rcpp::as<std::string>(names[i]))) {
            NumericVector param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] =
                as<arma::colvec>(param);
          } else if (isNumericMatrixInList(recipe_from_R(list_params),
                                           Rcpp::as<std::string>(names[i]))) {
            arma::mat param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          } else {
            stop("Initial value must be a float/double, vector or matrix.");
          }
        }
      } else {
        Rcpp::stop("List of initial values has empty entries.");
      }

      output.push_back(current_parameters);
    } else {
      Rcpp::stop("List of initial values must be a list (over chains) of "
                 "lists, each of which contain the initial parameters.");
    }
  }

  return output;
}

Parameters make_algorithm_parameters(const Rcpp::List &algorithm_parameter_list) {

  return list_to_parameters(algorithm_parameter_list);
}


ImportanceSampler *get_importance_sampler(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_importance_points, bool parallel_in, size_t grain_size_in,
    const std::string &results_name_in, size_t *seed,
    std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {
  // May need to alter for cases where the likelihood needs to be tuned
  // automatically (e.g. in ABC).

  // Check if the prior is the proposal: affects what llhd_estimators we
  // include.

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  IndependentProposalKernel *proposal_in;
  bool proposal_is_evaluated_in;

  RecipeValue sequencer_info = get_smc_sequencer_info(model, parameters, RecipeValue::list());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;

  if (model.contains("is_proposal")) {

    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, true,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = true;
  } else {
    inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, false,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = false;
  }


  return new ImportanceSampler(
      rng, seed, the_data, algorithm_parameters, number_of_importance_points,
      "", likelihood_estimators, proposal_in, without_cancelled_index, proposal_is_evaluated_in, true,
      true, true, false, parallel_in, grain_size_in, "");
}

// Output options set by the R wrappers as elements of the recipe.
static void apply_output_options(const RecipeValue &model, SMC *alg) {
  if (model.contains("write_factor_output"))
    alg->set_write_factor_output((model["write_factor_output"]).as_bool());
  // testing switch: never keep a generation in flat form (dev/plan.md 2.6a)
  if (model.contains("force_rich_generations"))
    alg->set_force_rich_generations((model["force_rich_generations"]).as_bool());
  // threads for parallel runs (set by the R wrappers from RcppParallel)
  if (model.contains("number_of_threads"))
    alg->set_number_of_threads(size_t((model["number_of_threads"]).as_double()));
  // progress bar for this run (not for inner estimators)
  if (model.contains("show_progress"))
    alg->set_show_progress((model["show_progress"]).as_bool(),
                           model.contains("progress_counts_iterations") &&
                           (model["progress_counts_iterations"]).as_bool());
  // waste-free SMC (SMC_with_MCMC only)
  if (model.contains("waste_free_chains"))
    alg->set_waste_free_chains(size_t((model["waste_free_chains"]).as_double()));
}

// After a run with parallel = TRUE (dev/plan.md 2.6a, G4).
static void warn_if_parallel_refused(const SMC *alg) {
  if (alg->parallel_was_refused())
    Rcpp::warning("parallel = TRUE: part of this run was done sequentially, because "
                  "it was not on the chunked path or a factor estimates its likelihood "
                  "with an inner algorithm. The results are the same as a sequential run.");
}

// [[Rcpp::export]]
double do_importance_sampler(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                             const Rcpp::List &algorithm_parameter_list_R,
                             size_t number_of_importance_points,
                             bool parallel_in, size_t grain_size_in,
                             const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_measurement_covariance_estimators;

  ImportanceSampler *alg = get_importance_sampler(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_importance_points, parallel_in, grain_size_in, results_name_in,
      &seed, data_created_in_get_likelihood_estimators,
      data_created_in_measurement_covariance_estimators);

  apply_output_options(model, alg);
  // nothing will be written: the factor outputs need not be kept
  if (strcmp(results_name_in.c_str(), "") == 0)
    alg->set_write_factor_output(false);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run();
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (strcmp(results_name_in.c_str(), "") != 0)
    output->write(results_name_in.c_str());

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

// [[Rcpp::export]]
double do_importance_sampler_with_fixed_params(
    const Rcpp::List &model_R, const Rcpp::List &parameters_R,
    const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &fixed_parameter_list,
    size_t number_of_importance_points, bool parallel_in, size_t grain_size_in,
    const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_measurement_covariance_estimators;

  ImportanceSampler *alg = get_importance_sampler(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_importance_points, parallel_in, grain_size_in, results_name_in,
      &seed, data_created_in_get_likelihood_estimators,
      data_created_in_measurement_covariance_estimators);

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run(list_to_parameters(fixed_parameter_list));
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (strcmp(results_name_in.c_str(), "") != 0)
    output->write(results_name_in.c_str());

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

// [[Rcpp::export]]
void do_mcmc(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
             const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &initial_values,
             const Rcpp::List &mcmc_termination_method_R,
             const Rcpp::List &mcmc_weights_method_R, size_t number_of_chains,
             bool parallel_in, size_t grain_size_in,
             const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);

  // could specify with mhproposalkernels
  // if more than one, need to specify type of sweep
  // for each need to specify if g, m or mh proposal, also which factors
  // involved

  // either use initial values list if supplied, else use prior

  // ilike:: includes any of the options we have made

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  // May need to alter for cases where the likelihood needs to be tuned
  // automatically (e.g. in ABC).

  // Check if the prior is the proposal: affects what llhd_estimators we
  // include.

  RecipeValue sequencer_info = get_smc_sequencer_info(model, parameters, RecipeValue::list());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, true,
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, any_annealing,
      factors_affected_by_smc_sequence,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  MCMC *the_mcmc =
      make_mcmc(model, parameters, &the_data, mcmc_termination_method,
                mcmc_weights_method, full_index);

  SMCMCMCMove *alg;

  if (initial_values.size() == 0) {

    // no initial values: the chains start at independent draws from the
    // recipe's is_proposal if it has one, otherwise from the prior
    IndependentProposalKernel *proposal_in =
        model.contains("is_proposal")
            ? get_proposal(model, parameters, &the_data)
            : get_prior_as_proposal(model, parameters);

    alg = new SMCMCMCMove(
        &rng, &seed, &the_data, algorithm_parameters, number_of_chains, 2, 2,
        the_mcmc, likelihood_estimators, proposal_in, without_cancelled_index,
        full_index, true, parallel_in, grain_size_in, "");
  } else {
    std::vector<Parameters> initial_points =
        make_initial_points(initial_values);

    arma::colvec log_probabilities_of_initial_values(initial_points.size());

    alg = new SMCMCMCMove(&rng, &seed, &the_data, algorithm_parameters, 2, 2,
                          the_mcmc, likelihood_estimators, initial_points,
                          log_probabilities_of_initial_values,
                          without_cancelled_index, full_index, true,
                          parallel_in, grain_size_in, "");
  }

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run();
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (strcmp(results_name_in.c_str(), "") != 0)
    output->write(results_name_in.c_str());

  delete full_index; // MCMC moves keep their own copies of it
}

// [[Rcpp::export]]
void do_mcmc_with_fixed_params(
    const Rcpp::List &model_R, const Rcpp::List &parameters_R,
    const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &fixed_parameter_list,
    const Rcpp::List &initial_values, const Rcpp::List &mcmc_termination_method_R,
    const Rcpp::List &mcmc_weights_method_R, size_t number_of_chains, bool parallel_in,
    size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);

  // could specify with mhproposalkernels
  // if more than one, need to specify type of sweep
  // for each need to specify if g, m or mh proposal, also which factors
  // involved

  // either use initial values list if supplied, else use prior

  // ilike:: includes any of the options we have made

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  // May need to alter for cases where the likelihood needs to be tuned
  // automatically (e.g. in ABC).

  // Check if the prior is the proposal: affects what llhd_estimators we
  // include.

  RecipeValue sequencer_info = get_smc_sequencer_info(model, parameters, RecipeValue::list());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, true,
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, any_annealing,
      factors_affected_by_smc_sequence,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  MCMC *the_mcmc =
      make_mcmc(model, parameters, &the_data, mcmc_termination_method,
                mcmc_weights_method, full_index);

  SMCMCMCMove *alg;

  if (initial_values.size() == 0) {

    // no initial values: the chains start at independent draws from the
    // recipe's is_proposal if it has one, otherwise from the prior
    IndependentProposalKernel *proposal_in =
        model.contains("is_proposal")
            ? get_proposal(model, parameters, &the_data)
            : get_prior_as_proposal(model, parameters);

    alg = new SMCMCMCMove(
        &rng, &seed, &the_data, algorithm_parameters, number_of_chains, 2, 2,
        the_mcmc, likelihood_estimators, proposal_in, without_cancelled_index,
        full_index, true, parallel_in, grain_size_in, "");
  } else {
    std::vector<Parameters> initial_points =
        make_initial_points(initial_values);

    arma::colvec log_probabilities_of_initial_values(initial_points.size());

    alg = new SMCMCMCMove(&rng, &seed, &the_data, algorithm_parameters, 2, 2,
                          the_mcmc, likelihood_estimators, initial_points,
                          log_probabilities_of_initial_values,
                          without_cancelled_index, full_index, true,
                          parallel_in, grain_size_in, "");
  }

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run(list_to_parameters(fixed_parameter_list));
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (strcmp(results_name_in.c_str(), "") != 0)
    output->write(results_name_in.c_str());

  delete full_index; // MCMC moves keep their own copies of it
}

SMCMCMCMove *get_smc_mcmc_move(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_particles, const RecipeValue &mcmc_termination_method,
    const RecipeValue &mcmc_weights_method, const RecipeValue &smc_sequencer_method,
    const RecipeValue &adaptive_target_method, const RecipeValue &smc_termination_method,
    size_t smc_iterations_to_store, bool write_to_file_at_each_iteration,
    bool parallel_in, size_t grain_size_in, const std::string &results_name_in,
    size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {
  std::vector<LikelihoodEstimator *> likelihood_estimators;
  IndependentProposalKernel *proposal_in;
  bool proposal_is_evaluated_in;

  RecipeValue sequencer_info =
      get_smc_sequencer_info(model, parameters, smc_sequencer_method);
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;

  bool any_annealing = false;

  if (model.contains("is_proposal")) {
    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, true,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = true;
  } else {
    inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, false,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = false;
  }

  // note that this will always result in the temperature being adapted in the
  // SMC algorithm (not smcfixed)


  MCMC *the_mcmc =
      make_mcmc(model, parameters, the_data, mcmc_termination_method,
                mcmc_weights_method, full_index);

  SMCCriterion *resampling_criterion =
      get_resampling_method(model, parameters, number_of_particles);

  SMCCriterion *adaptive_target_criterion = get_adaptive_target_method(
      model, parameters, adaptive_target_method, number_of_particles);

  SMCTermination *smc_termination =
      get_smc_termination_method(model, parameters, smc_termination_method);

  SMCMCMCMove *alg;
  if (write_to_file_at_each_iteration) {
    alg = new SMCMCMCMove(
        rng, seed, the_data, algorithm_parameters, number_of_particles,
        smc_iterations_to_store, smc_iterations_to_store, the_mcmc,
        resampling_criterion, adaptive_target_criterion, sequencer_info[3],
        smc_termination, sequencer_variables, sequencer_schedules,
        likelihood_estimators, proposal_in, without_cancelled_index, full_index,
        proposal_is_evaluated_in, true, true, false, false, parallel_in,
        grain_size_in, results_name_in.c_str());
  } else {
    alg = new SMCMCMCMove(
        rng, seed, the_data, algorithm_parameters, number_of_particles,
        smc_iterations_to_store, smc_iterations_to_store, the_mcmc,
        resampling_criterion, adaptive_target_criterion, sequencer_info[3],
        smc_termination, sequencer_variables, sequencer_schedules,
        likelihood_estimators, proposal_in, without_cancelled_index, full_index,
        proposal_is_evaluated_in, true, true, false, false, parallel_in,
        grain_size_in, "");
  }
  delete full_index; // MCMC moves keep their own copies of it
  return alg;
}

// [[Rcpp::export]]
double do_smc_mcmc_move(
    const Rcpp::List &model_R, const Rcpp::List &parameters_R,
    const Rcpp::List &algorithm_parameter_list_R, size_t number_of_particles,
    const Rcpp::List &mcmc_termination_method_R, const Rcpp::List &mcmc_weights_method_R,
    const Rcpp::List &smc_sequencer_method_R, const Rcpp::List &adaptive_target_method_R,
    const Rcpp::List &smc_termination_method_R, size_t smc_iterations_to_store,
    bool write_to_file_at_each_iteration, bool parallel_in,
    size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);
  RecipeValue smc_sequencer_method = recipe_from_R(smc_sequencer_method_R);
  RecipeValue adaptive_target_method = recipe_from_R(adaptive_target_method_R);
  RecipeValue smc_termination_method = recipe_from_R(smc_termination_method_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  SMCMCMCMove *alg = get_smc_mcmc_move(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_particles, mcmc_termination_method, mcmc_weights_method,
      smc_sequencer_method, adaptive_target_method, smc_termination_method,
      smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in,
      grain_size_in, results_name_in, &seed,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run();
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

// [[Rcpp::export]]
double do_smc_mcmc_move_with_fixed_params(
                        const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                        const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &fixed_parameter_list, size_t number_of_particles,
                        const Rcpp::List &mcmc_termination_method_R, const Rcpp::List &mcmc_weights_method_R,
                        const Rcpp::List &smc_sequencer_method_R, const Rcpp::List &adaptive_target_method_R,
                        const Rcpp::List &smc_termination_method_R, size_t smc_iterations_to_store,
                        bool write_to_file_at_each_iteration, bool parallel_in,
                        size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);
  RecipeValue smc_sequencer_method = recipe_from_R(smc_sequencer_method_R);
  RecipeValue adaptive_target_method = recipe_from_R(adaptive_target_method_R);
  RecipeValue smc_termination_method = recipe_from_R(smc_termination_method_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  SMCMCMCMove *alg = get_smc_mcmc_move(
                                       &rng, &the_data, model, parameters, algorithm_parameters,
                                       number_of_particles, mcmc_termination_method, mcmc_weights_method,
                                       smc_sequencer_method, adaptive_target_method, smc_termination_method,
                                       smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in,
                                       grain_size_in, results_name_in, &seed,
                                       data_created_in_get_likelihood_estimators,
                                       data_created_in_get_measurement_covariance_estimators);

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run(list_to_parameters(fixed_parameter_list));
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

EnsembleKalmanInversion *get_enki(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_ensemble_members, const RecipeValue &mcmc_termination_method,
    const RecipeValue &mcmc_weights_method, const RecipeValue &enk_sequencer_method,
    const RecipeValue &adaptive_target_method, const RecipeValue &enk_termination_method,
    size_t enk_iterations_to_store, bool write_to_file_at_each_iteration,
    bool parallel_in, size_t grain_size_in, const std::string &results_name_in,
    size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {
  std::vector<LikelihoodEstimator *> likelihood_estimators;
  std::vector<MeasurementCovarianceEstimator *> estimators;
  IndependentProposalKernel *proposal_in;

  RecipeValue sequencer_info =
      get_smc_sequencer_info(model, parameters, enk_sequencer_method,
                             BuiltinRole::EnKSequence);
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;

  bool any_annealing = false;

  proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

  std::shared_ptr<Transform> transform = NULL;

  std::vector<size_t> enk_indices =
      get_enk_likelihood_indices(model, parameters);

  likelihood_estimators = get_likelihood_estimators(
      rng, seed, the_data, model, parameters, algorithm_parameters, true,
      sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
      without_cancelled_index, without_priors_index, full_index, any_annealing,
      factors_affected_by_smc_sequence,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  estimators = get_measurement_covariance_estimators(
      rng, seed, the_data, model, parameters, enk_indices, transform,
      data_created_in_get_measurement_covariance_estimators);

  // note that this will always result in the temperature being adapted in the
  // SMC algorithm (not smcfixed)

  MCMC *the_mcmc =
      make_mcmc(model, parameters, the_data, mcmc_termination_method,
                mcmc_weights_method, full_index);

  SMCCriterion *adaptive_target_criterion = get_adaptive_target_method(
      model, parameters, adaptive_target_method, number_of_ensemble_members);

  SMCTermination *enk_termination =
      get_smc_termination_method(model, parameters, enk_termination_method);

  EnsembleShifter *shifter = get_enk_shifter_method(model, parameters);

  EnsembleKalmanInversion *alg;
  if (write_to_file_at_each_iteration) {
    alg = new EnsembleKalmanInversion(
        rng, seed, the_data, number_of_ensemble_members,
        enk_iterations_to_store, shifter, adaptive_target_criterion,
        sequencer_info[3], enk_termination,
        sequencer_variables[sequencer_variables.size() - 1],
        sequencer_schedules[sequencer_schedules.size() - 1], proposal_in,
        likelihood_estimators, estimators, transform, 1.0, 0, parallel_in,
        grain_size_in, results_name_in.c_str());
  } else {
    alg = new EnsembleKalmanInversion(
        rng, seed, the_data, number_of_ensemble_members,
        enk_iterations_to_store, shifter, adaptive_target_criterion,
        sequencer_info[3], enk_termination,
        sequencer_variables[sequencer_variables.size() - 1],
        sequencer_schedules[sequencer_schedules.size() - 1], proposal_in,
        likelihood_estimators, estimators, transform, 1.0, 0, parallel_in,
        grain_size_in, "");
  }
  delete full_index; // MCMC moves keep their own copies of it
  return alg;
}

// [[Rcpp::export]]
double
do_enki(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
        const Rcpp::List &algorithm_parameter_list_R, size_t number_of_ensemble_members,
        const Rcpp::List &mcmc_termination_method_R, const Rcpp::List &mcmc_weights_method_R,
        const Rcpp::List &enk_sequencer_method_R, const Rcpp::List &adaptive_target_method_R,
        const Rcpp::List &enk_termination_method_R, size_t enk_iterations_to_store,
        bool write_to_file_at_each_iteration, bool parallel_in,
        size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);
  RecipeValue enk_sequencer_method = recipe_from_R(enk_sequencer_method_R);
  RecipeValue adaptive_target_method = recipe_from_R(adaptive_target_method_R);
  RecipeValue enk_termination_method = recipe_from_R(enk_termination_method_R);

  RandomNumberGenerator rng;

  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  EnsembleKalmanInversion *alg = get_enki(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_ensemble_members, mcmc_termination_method, mcmc_weights_method,
      enk_sequencer_method, adaptive_target_method, enk_termination_method,
      enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in,
      grain_size_in, results_name_in, &seed,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  EnsembleKalmanOutput *output = alg->run();

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  delete output;
  delete alg;

  return log_likelihood;
}

// MFDS: the ensemble Kalman method of Huang, Huang, Reich & Stuart
// (arXiv-2204.04386v2; dev/plan.md 3.16). The posterior is the fixed point of
// a filtering problem: each iteration inflates the ensemble about its mean
// (prediction) and updates it with the augmented observation [y; prior mean]
// of [simulated data; theta], whose noise is scaled by 1/delta_t (analysis).
// The prior must be built-in Gaussian; the likelihood is a simulator
// (simulate_data_model, data_variable, optional summary_statistics). The
// shifter (enk_shifter) is stochastic (EKI, default) or adjustment (EAKI).
// There is no evidence estimate.
// [[Rcpp::export]]
void do_enkmfds(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                size_t number_of_ensemble_members, size_t number_of_iterations,
                double delta_t, size_t iterations_to_store,
                const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);

  if (!((delta_t > 0.0) && (delta_t < 1.0)))
    stop("MFDS: delta_t must be between 0 and 1.");
  if (number_of_iterations < 1)
    stop("MFDS: number_of_iterations must be at least 1.");
  if (number_of_ensemble_members < 2)
    stop("MFDS: at least two ensemble members are needed.");

  RandomNumberGenerator rng;
  Data the_data = get_data(model);

  // the Gaussian prior: its means are the artificial observation of theta
  Parameters prior_means;
  Parameters prior_covariances;
  get_gaussian_prior_moments(model, parameters, prior_means, prior_covariances);

  // the simulator of the data
  SimulateModelPtr simulate_data_model = NULL;
  std::vector<std::string> data_variables;
  std::shared_ptr<Transform> summary_statistics = NULL;
  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];
    for (size_t i = 0; i < factors.size(); ++i) {
      RecipeValue current_factor = factors[i];
      if (!current_factor.contains("simulate_data_model"))
        continue;
      if (simulate_data_model != NULL)
        stop("MFDS: only one simulate_data_model is allowed.");
      simulate_data_model =
          current_factor["simulate_data_model"].function_pointer<SimulateModelPtr>();
      if (!current_factor.contains("data_variable"))
        stop("MFDS: the factor with simulate_data_model must also have data_variable.");
      data_variables = split(current_factor["data_variable"].function_pointer<StringPtr>()(), ';');
      if (current_factor.contains("summary_statistics"))
        summary_statistics = std::make_shared<Transform>(
            current_factor["summary_statistics"].function_pointer<SummaryStatisticsPtr>());
    }
  }
  if (simulate_data_model == NULL)
    stop("MFDS: the model needs a simulate_data_model (and data_variable) for the likelihood.");

  // the initial ensemble, drawn from the prior (with the streams' convention
  // of the ensemble workers: seed, then advance the seed by the ensemble size)
  std::unique_ptr<IndependentProposalKernel> prior(get_prior_as_simulate_only_proposal(model, parameters));
  rng.seed(seed);
  std::vector<Parameters> initial_points;
  initial_points.reserve(number_of_ensemble_members);
  for (size_t i = 0; i < number_of_ensemble_members; ++i)
    initial_points.push_back(prior->independent_simulate(rng));
  seed = seed + number_of_ensemble_members;

  EnsembleShifter *shifter = get_enk_shifter_method(model, parameters);
  if (dynamic_cast<SquareRootEnsembleShifter*>(shifter) != NULL) {
    delete shifter;
    stop("MFDS: enk_shifter ilike::sqrt is not available for MFDS; use ilike::stochastic or ilike::adjustment.");
  }

  std::unique_ptr<EnsembleKalmanMFDS> alg(new EnsembleKalmanMFDS(
      &rng, &seed, &the_data, iterations_to_store, initial_points, shifter, delta_t,
      number_of_iterations, prior_means, prior_covariances, simulate_data_model,
      2, // the joint update of the augmented observation
      summary_statistics, NULL, data_variables, false, 100000, results_name_in));

  std::unique_ptr<EnsembleKalmanOutput> output(alg->run());
}

KalmanFilter *get_kalman_filter(Data *the_data, const RecipeValue &model,
                                const RecipeValue &parameters,
                                const Parameters &algorithm_parameters,
                                size_t kf_iterations_to_store,
                                bool write_to_file_at_each_iteration,
                                const std::string &results_name_in) {
  RecipeValue prior_mean_and_covariance =
      get_prior_mean_and_covariance(model, parameters);

  arma::colvec prior_mean_in = prior_mean_and_covariance[0];
  arma::mat prior_covariance_in = prior_mean_and_covariance[1];

  KalmanPredictor *predictor_in = get_kalman_predictor(model, parameters);

  std::vector<KalmanUpdater *> updaters_in =
      get_kalman_updaters(model, parameters);

  RecipeValue filtering_info = get_ssm_info(model, parameters);

  std::string index_name_in = filtering_info[0];
  size_t first_index_in = filtering_info[1];
  size_t last_index_in = filtering_info[2];
  std::string time_name_in = filtering_info[3];
  double initial_time_in = filtering_info[4];
  std::string time_diff_name_in = filtering_info[5];
  double update_time_step_in = filtering_info[6];
  size_t predictions_per_update_in = filtering_info[7];

  if (first_index_in > last_index_in)
    stop("Cannot construct a filter where the first index is bigger than the "
         "last.");

  for (auto i = updaters_in.begin(); i != updaters_in.end(); ++i) {
    if (((*i)->get_measurement_variable() == "") &&
        (last_index_in >= the_data->min_n_rows()))
      stop("Last index goes past the end of the data.");

    if (((*i)->get_measurement_variable() != "") &&
        (last_index_in >= (*the_data)[(*i)->get_measurement_variable()].n_rows))
      stop("Last index goes past the end of the data.");
  }

  // KalmanPredictor* predictor_in = new
  // ExactKalmanPredictor(transition_model_A,
  //                                                          transition_model_Q);
  // KalmanUpdater* updater_in = new ExactKalmanUpdater(measurement_model_H,
  //                                                    measurement_model_R);


  if (write_to_file_at_each_iteration) {
    return new KalmanFilter(the_data, kf_iterations_to_store,
                            // state_name_in,
                            prior_mean_in, prior_covariance_in, index_name_in,
                            time_name_in, time_diff_name_in,
                            // measurement_names_in,
                            first_index_in, last_index_in,
                            predictions_per_update_in, update_time_step_in,
                            initial_time_in, true, predictor_in, updaters_in,
                            true, results_name_in);
  } else {
    return new KalmanFilter(
        the_data, kf_iterations_to_store, prior_mean_in, prior_covariance_in,
        index_name_in, time_name_in, time_diff_name_in, first_index_in,
        last_index_in, predictions_per_update_in, update_time_step_in,
        initial_time_in, true, predictor_in, updaters_in, true, "");
  }
}

// [[Rcpp::export]]
double do_kalman_filter(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                        const Rcpp::List &algorithm_parameter_list_R,
                        size_t kf_iterations_to_store,
                        bool write_to_file_at_each_iteration,
                        const std::string &results_name_in) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  Data the_data = get_data(model);

  KalmanFilter *alg =
      get_kalman_filter(&the_data, model, parameters, algorithm_parameters,
                        kf_iterations_to_store, write_to_file_at_each_iteration,
                        results_name_in.c_str());

  KalmanFilterOutput *output = alg->run();

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  delete output;
  delete alg;

  return log_likelihood;
}

EnsembleKalmanFilter *get_ensemble_kalman_filter(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_ensemble_members, EnsembleShifter *shifter,
    // const RecipeValue &enk_likelihood_index_method,
    // const RecipeValue &enk_shifter_method,
    size_t enk_iterations_to_store, bool write_to_file_at_each_iteration,
    bool parallel_in, size_t grain_size_in, const std::string &results_name_in,
    size_t *seed,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {

  std::vector<MeasurementCovarianceEstimator *> estimators;

  IndependentProposalKernel *proposal_in;
  proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

  std::vector<size_t> enk_likelihood_indices =
      get_enk_likelihood_indices(model, parameters);

  bool smcfixed_flag = true;

  RecipeValue filtering_info = get_ssm_info(model, parameters);

  std::string index_name_in = filtering_info[0];
  size_t first_index_in = filtering_info[1];
  size_t last_index_in = filtering_info[2];
  std::string time_name_in = filtering_info[3];
  double initial_time_in = filtering_info[4];
  std::string time_diff_name_in = filtering_info[5];
  double update_time_step_in = filtering_info[6];
  size_t predictions_per_update_in = filtering_info[7];

  if (first_index_in > last_index_in)
    stop("Cannot construct a filter where the first index is bigger than the "
         "last.");

  ProposalKernel *transition = get_transition_model(model, parameters);

  std::vector<MeasurementCovarianceEstimator *>
      measurement_covariance_estimators_in;

  std::shared_ptr<Transform> transform = NULL;

  std::vector<size_t> enk_indices =
      get_enk_likelihood_indices(model, parameters);

  measurement_covariance_estimators_in = get_measurement_covariance_estimators(
      rng, seed, the_data, model, parameters, enk_indices, transform,
      data_created_in_get_measurement_covariance_estimators);

  for (auto i = measurement_covariance_estimators_in.begin();
       i != measurement_covariance_estimators_in.end(); ++i) {
    std::vector<std::string> measurement_variables =
        (*i)->get_measurement_variables();

    for (auto j = measurement_variables.begin();
         j != measurement_variables.end(); ++j) {
      if ((*j == "") && (last_index_in >= the_data->min_n_rows()))
        stop("Last index goes past the end of the data.");

      if ((*j != "") && (last_index_in >= (*the_data)[*j].n_rows))
        stop("Last index goes past the end of the data.");
    }
  }

  if (write_to_file_at_each_iteration) {
    return new EnsembleKalmanFilter(
        rng, seed, the_data, enk_iterations_to_store, index_name_in,
        time_name_in, time_diff_name_in, first_index_in, last_index_in,
        predictions_per_update_in, update_time_step_in, initial_time_in,
        number_of_ensemble_members, shifter, NULL, proposal_in, transition,
        measurement_covariance_estimators_in, smcfixed_flag, true, parallel_in,
        grain_size_in, results_name_in.c_str());
  } else {
    return new EnsembleKalmanFilter(
        rng, seed, the_data, enk_iterations_to_store, index_name_in,
        time_name_in, time_diff_name_in, first_index_in, last_index_in,
        predictions_per_update_in, update_time_step_in, initial_time_in,
        number_of_ensemble_members, shifter, NULL, proposal_in, transition,
        measurement_covariance_estimators_in, smcfixed_flag, true, parallel_in,
        grain_size_in, "");
  }
}

// [[Rcpp::export]]
double do_ensemble_kalman_filter(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                                 const Rcpp::List &algorithm_parameter_list_R,
                                 size_t number_of_ensemble_members,
                                 size_t enk_iterations_to_store,
                                 bool write_to_file_at_each_iteration,
                                 bool parallel_in, size_t grain_size_in,
                                 const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);

  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  EnsembleShifter *shifter = get_enk_shifter_method(model, parameters);

  EnsembleKalmanFilter *alg = get_ensemble_kalman_filter(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_ensemble_members, shifter, enk_iterations_to_store,
      write_to_file_at_each_iteration, parallel_in, grain_size_in,
      results_name_in, &seed,
      data_created_in_get_measurement_covariance_estimators);

  EnsembleKalmanOutput *output = alg->run();

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  delete output;
  delete alg;

  return log_likelihood;
}

// [[Rcpp::export]]
double do_ensemble_kalman_filter_with_fixed_params(
    const Rcpp::List &model_R, const Rcpp::List &parameters_R,
    const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &fixed_parameter_list,
    size_t number_of_ensemble_members, size_t enk_iterations_to_store,
    bool write_to_file_at_each_iteration, bool parallel_in,
    size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);

  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  EnsembleShifter *shifter = get_enk_shifter_method(model, parameters);

  EnsembleKalmanFilter *alg = get_ensemble_kalman_filter(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_ensemble_members, shifter, enk_iterations_to_store,
      write_to_file_at_each_iteration, parallel_in, grain_size_in,
      results_name_in, &seed,
      data_created_in_get_measurement_covariance_estimators);

  EnsembleKalmanOutput *output =
      alg->run(list_to_parameters(fixed_parameter_list));

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  delete output;
  delete alg;

  return log_likelihood;
}

ParticleFilter *get_particle_filter(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_particles, size_t smc_iterations_to_store,
    bool write_to_file_at_each_iteration, bool parallel_in,
    size_t grain_size_in, const std::string &results_name_in, size_t *seed,
    std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {
  std::vector<LikelihoodEstimator *> likelihood_estimators;
  IndependentProposalKernel *proposal_in;
  bool proposal_is_evaluated_in;

  VectorIndex *evaluated_in_initial_weight_update = NULL;
  VectorIndex *evaluated_in_pf_weight_update = NULL;
  VectorIndex *full_index = NULL;

  std::vector<std::string> sequencer_types;
  std::vector<std::string> sequencer_variables;
  std::vector<std::vector<double>> sequencer_schedules;

  std::vector<int> factors_affected_by_smc_sequence;

  bool any_annealing;

  if (model.contains("is_proposal")) {

    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, true,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        full_index, any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = true;
  } else {

    inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, false,
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        full_index, any_annealing, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = false;
  }

  bool smcfixed_flag = true;

  SMCCriterion *resampling_criterion =
      get_resampling_method(model, parameters, number_of_particles);

  RecipeValue filtering_info = get_ssm_info(model, parameters);

  std::string index_name_in = filtering_info[0];
  size_t first_index_in = filtering_info[1];
  size_t last_index_in = filtering_info[2];
  std::string time_name_in = filtering_info[3];
  double initial_time_in = filtering_info[4];
  std::string time_diff_name_in = filtering_info[5];
  double update_time_step_in = filtering_info[6];
  size_t predictions_per_update_in = filtering_info[7];

  if (first_index_in > last_index_in)
    stop("Cannot construct a filter where the first index is bigger than the "
         "last.");

  // Commented out with the intention to replace - issue is that we no longer
  // store the measurement name.

  ProposalKernel *transition_model = get_transition_model(model, parameters);

  ProposalKernel *transition_proposal =
      get_transition_proposal(model, parameters);

  bool transition_proposal_is_evaluated_in;
  if (transition_proposal == NULL) {
    inform_once("No transition proposal specified; using transition model.");

    transition_proposal_is_evaluated_in = false;

    transition_proposal = transition_model->proposal_kernel_duplicate();
  } else {
    transition_proposal_is_evaluated_in = true;

    if (!transition_model->can_be_evaluated()) {
      stop("No method supplied for the evaluation of the transition model.");
    }
  }

  if (write_to_file_at_each_iteration) {
    return new ParticleFilter(
        rng, seed, the_data, Parameters(), number_of_particles,
        smc_iterations_to_store, smc_iterations_to_store, index_name_in,
        time_name_in, time_diff_name_in, first_index_in, last_index_in,
        predictions_per_update_in, update_time_step_in, initial_time_in,
        resampling_criterion, likelihood_estimators, proposal_in,
        transition_model, transition_proposal,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        proposal_is_evaluated_in, transition_proposal_is_evaluated_in,
        smcfixed_flag, true, false, parallel_in, grain_size_in,
        results_name_in.c_str());
  } else {
    return new ParticleFilter(
        rng, seed, the_data, Parameters(), number_of_particles,
        smc_iterations_to_store, smc_iterations_to_store, index_name_in,
        time_name_in, time_diff_name_in, first_index_in, last_index_in,
        predictions_per_update_in, update_time_step_in, initial_time_in,
        resampling_criterion, likelihood_estimators, proposal_in,
        transition_model, transition_proposal,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        proposal_is_evaluated_in, transition_proposal_is_evaluated_in,
        smcfixed_flag, true, false, parallel_in, grain_size_in, "");
  }
}

// [[Rcpp::export]]
double do_particle_filter(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                          const Rcpp::List &algorithm_parameter_list_R,
                          size_t number_of_particles,
                          size_t smc_iterations_to_store,
                          bool write_to_file_at_each_iteration,
                          bool parallel_in, size_t grain_size_in,
                          const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);

  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  ParticleFilter *alg = get_particle_filter(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_particles, smc_iterations_to_store,
      write_to_file_at_each_iteration, parallel_in, grain_size_in,
      results_name_in, &seed, data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run();
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

// [[Rcpp::export]]
double do_particle_filter_with_fixed_params(
    const Rcpp::List &model_R, const Rcpp::List &parameters_R,
    const Rcpp::List &algorithm_parameter_list_R, const Rcpp::List &fixed_parameter_list,
    size_t number_of_particles, size_t smc_iterations_to_store,
    bool write_to_file_at_each_iteration, bool parallel_in,
    size_t grain_size_in, const std::string &results_name_in, size_t seed) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  RandomNumberGenerator rng;

  Data the_data = get_data(model);

  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  ParticleFilter *alg = get_particle_filter(
      &rng, &the_data, model, parameters, algorithm_parameters,
      number_of_particles, smc_iterations_to_store,
      write_to_file_at_each_iteration, parallel_in, grain_size_in,
      results_name_in, &seed, data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  apply_output_options(model, alg);

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run(list_to_parameters(fixed_parameter_list));
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);

  if (!write_to_file_at_each_iteration) {
    if (strcmp(results_name_in.c_str(), "") != 0)
      output->write(results_name_in.c_str());
  }

  double log_likelihood = output->log_likelihood;

  return log_likelihood;
}

// Factor numbers given by the user in evaluate_log are 1-based (as in .ilike
// files); C++ indexes the factors from 0.
VectorIndex* user_factor_index(const IntegerVector &index_in,
                               size_t number_of_factors)
{
  std::vector<size_t> indices_in;
  for (R_xlen_t k = 0; k < index_in.size(); ++k) {
    if ((index_in[k] == NA_INTEGER) || (index_in[k] < 1) ||
        (size_t(index_in[k]) > number_of_factors))
      Rcpp::stop("evaluate_log: index must contain factor numbers between 1 and " +
                 std::to_string(number_of_factors) + ".");
    indices_in.push_back(size_t(index_in[k] - 1));
  }
  return new VectorIndex(indices_in);
}

// [[Rcpp::export]]
double do_evaluate_log(const Rcpp::List &model_R,
                       const Rcpp::List &model_parameters_R,
                       const Rcpp::List &algorithm_parameter_list_R,
                       const Rcpp::List &parameters_in,
                       const IntegerVector &index_in,
                       const std::string &results_name_in,
                       size_t seed)
{
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue model_parameters = recipe_from_R(model_parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);


  Parameters parameters = list_to_parameters(parameters_in);

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  std::vector<std::string> sequencer_types;
  std::vector<std::string> sequencer_variables;
  std::vector<std::vector<double>> sequencer_schedules;

  std::vector<int> factors_affected_by_smc_sequence;

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, model_parameters, algorithm_parameters, true,
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, any_annealing,
      factors_affected_by_smc_sequence,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  Factors* factors = new VectorFactors(likelihood_estimators);

  Index* index;
  if (index_in.size() == 0) {
    index = full_index;
  }
  else
  {
    index = user_factor_index(index_in, likelihood_estimators.size());
  }

  factors->setup(parameters);

  // augment with conditioned on...
  //parameters.merge_with_fixed(conditioned_on_parameters);

  std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

  const std::vector<const ProposalKernel*>* proposals_to_transform_for_in = NULL;
  const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in = NULL;
  Particle particle(parameters,
                    factors,
                    proposals_to_transform_for_in,
                    proposals_to_find_gradient_for_in);

  particle.evaluate_smcfixed_part_of_likelihoods(index);
  double output = particle.evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);

  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;

  if (strcmp(results_name_in.c_str(), "") != 0)
  {
    std::string directory_name_in = results_name_in.c_str();
    std::string directory_name = directory_name_in + "_evaluate";

    if (!directory_exists(directory_name))
    {
      make_directory(directory_name);
    }

    auto h5f = h5_open_or_create(directory_name + "/output.h5");
    h5_append_double(*h5f, "log_likelihood", output);
    h5_append_double(*h5f, "time", elapsed_time.count());

    std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(1);

    particle.factor_variables->write_to_file(smc_iteration_directory,0);

  }

  return output;
}

// [[Rcpp::export]]
double do_evaluate_log_with_fixed_params(const Rcpp::List &model_R,
                       const Rcpp::List &model_parameters_R,
                       const Rcpp::List &algorithm_parameter_list_R,
                       const Rcpp::List &fixed_parameter_list,
                       const Rcpp::List &parameters_in,
                       const IntegerVector &index_in,
                       const std::string &results_name_in,
                       size_t seed)
{
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue model_parameters = recipe_from_R(model_parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);

  RandomNumberGenerator rng;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  std::vector<std::string> sequencer_types;
  std::vector<std::string> sequencer_variables;
  std::vector<std::vector<double>> sequencer_schedules;

  std::vector<int> factors_affected_by_smc_sequence;

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
                                                    &rng, &seed, &the_data, model, model_parameters, algorithm_parameters, true,
                                                    sequencer_types, sequencer_variables, sequencer_schedules, NULL,
                                                    without_cancelled_index, without_priors_index, full_index, any_annealing,
                                                    factors_affected_by_smc_sequence,
                                                    data_created_in_get_likelihood_estimators,
                                                    data_created_in_get_measurement_covariance_estimators);

  Factors* factors = new VectorFactors(likelihood_estimators);

  Index* index;
  if (index_in.size() == 0) {
    index = full_index;
  }
  else
  {
    index = user_factor_index(index_in, likelihood_estimators.size());
  }

  Parameters parameters = list_to_parameters(parameters_in);
  factors->setup(parameters);

  // augment with conditioned on...
  parameters.merge_with_fixed(list_to_parameters(fixed_parameter_list));

  std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();

  const std::vector<const ProposalKernel*>* proposals_to_transform_for_in = NULL;
  const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in = NULL;
  Particle particle(parameters,
                    factors,
                    proposals_to_transform_for_in,
                    proposals_to_find_gradient_for_in);

  particle.evaluate_smcfixed_part_of_likelihoods(index);
  double output = particle.evaluate_smcadaptive_part_given_smcfixed_likelihoods(index);

  std::chrono::high_resolution_clock::time_point end_time = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_time = end_time - start_time;

  if (strcmp(results_name_in.c_str(), "") != 0)
  {
    std::string directory_name_in = results_name_in.c_str();
    std::string directory_name = directory_name_in + "_evaluate";

    if (!directory_exists(directory_name))
    {
      make_directory(directory_name);
    }

    auto h5f = h5_open_or_create(directory_name + "/output.h5");
    h5_append_double(*h5f, "log_likelihood", output);
    h5_append_double(*h5f, "time", elapsed_time.count());

    std::string smc_iteration_directory = directory_name + "/iteration" + std::to_string(1);

    particle.factor_variables->write_to_file(smc_iteration_directory,0);

  }

  return output;
}

// Copy audit (dev/plan.md 3.5, CH1). Builds the likelihood estimators of a
// recipe as MCMC() does and, for each, evaluates it at a point simulated
// from the prior (the SMC sequence variables at the last value of their
// schedules), duplicates it, evaluates the copy, deletes the copy and
// evaluates the original again, all with the same seed: the three values
// must be identical. Then, if an MCMC move can be built, one move from the
// same state with the same seed by the move and by its copy (the copy then
// deleted, and the move run again). Returns a list of results; a step that
// stops is reported with its message. For tests (with AddressSanitizer).
static double audit_evaluate(LikelihoodEstimator *estimator,
                             const Parameters &parameters,
                             RandomNumberGenerator &rng, size_t *seed,
                             size_t seed_value) {
  *seed = seed_value;
  rng.seed(seed_value);
  LikelihoodEstimatorOutput *output = estimator->initialise(parameters);
  output->simulate(parameters);
  output->evaluate_smcfixed_part(parameters);
  output->evaluate_smcadaptive_part_given_smcfixed(parameters);
  double value = output->log_likelihood;
  delete output;
  return value;
}

// Internal check of shared_point_weights (used by the parameter estimators'
// fit from MoveOutputs, which R cannot reach): particle i's points are the
// rows of points_per_particle[[i]] (variable "x"). Returns the points (one
// per row) and their normalised log weights.
// [[Rcpp::export]]
Rcpp::List shared_point_weights_check(const Rcpp::List &points_per_particle,
                                const arma::colvec &normalised_log_weights) {

  std::vector< std::vector<Parameters> > points(points_per_particle.size());
  for (size_t i = 0; i < points.size(); ++i) {
    arma::mat m = Rcpp::as<arma::mat>(points_per_particle[i]);
    for (size_t r = 0; r < m.n_rows; ++r) {
      Parameters p;
      p["x"] = arma::colvec(m.row(r).t());
      points[i].push_back(p);
    }
  }
  std::vector<Parameters> all_points;
  arma::colvec all_normalised_log_weights;
  shared_point_weights(points, normalised_log_weights, all_points,
                       all_normalised_log_weights);
  return Rcpp::List::create(Named("points") = vector_of_parameters_to_mat("x", all_points),
                      Named("log_weights") = all_normalised_log_weights);
}

// [[Rcpp::export]]
Rcpp::List copy_audit(const Rcpp::List &model_R, const Rcpp::List &parameters_R,
                const Rcpp::List &algorithm_parameter_list_R,
                const Rcpp::List &mcmc_termination_method_R,
                const Rcpp::List &mcmc_weights_method_R, size_t seed_value) {
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue parameters = recipe_from_R(parameters_R);
  Parameters algorithm_parameters = make_algorithm_parameters(algorithm_parameter_list_R);
  RecipeValue mcmc_termination_method = recipe_from_R(mcmc_termination_method_R);
  RecipeValue mcmc_weights_method = recipe_from_R(mcmc_weights_method_R);

  RandomNumberGenerator rng;
  size_t seed = seed_value;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;

  RecipeValue sequencer_info = get_smc_sequencer_info(model, parameters, RecipeValue::list());
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];
  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() == 5) {
    std::vector<int> temp = sequencer_info[4];
    for (auto &k : temp) k = k - 1;
    factors_affected_by_smc_sequence = temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  bool any_annealing = false;
  std::vector<LikelihoodEstimator *> likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, true,
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, any_annealing,
      factors_affected_by_smc_sequence, data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);

  // a point from the prior, with the sequence variables at their final values
  IndependentProposalKernel *prior = get_prior_as_simulate_only_proposal(model, parameters);
  rng.seed(seed_value);
  Parameters point = prior->independent_simulate(rng);
  delete prior;
  for (size_t i = 0; i < sequencer_variables.size(); ++i) {
    if (sequencer_schedules[i].size() > 0)
      point[sequencer_variables[i]] = sequencer_schedules[i].back();
  }

  std::vector<std::string> messages;
  NumericVector original(likelihood_estimators.size(), NA_REAL);
  NumericVector copy(likelihood_estimators.size(), NA_REAL);
  NumericVector original_again(likelihood_estimators.size(), NA_REAL);
  for (size_t i = 0; i < likelihood_estimators.size(); ++i) {
    std::string message = "";
    try {
      // a first evaluation may set the estimator up (and use random numbers
      // doing so): compare from the second on
      audit_evaluate(likelihood_estimators[i], point, rng, &seed, seed_value);
      original[i] = audit_evaluate(likelihood_estimators[i], point, rng, &seed, seed_value);
      LikelihoodEstimator *duplicate = likelihood_estimators[i]->duplicate();
      copy[i] = audit_evaluate(duplicate, point, rng, &seed, seed_value);
      delete duplicate;
      original_again[i] = audit_evaluate(likelihood_estimators[i], point, rng, &seed, seed_value);
    } catch (std::exception &e) {
      message = e.what();
    }
    messages.push_back(message);
  }

  // one MCMC move by the move and by its copy
  std::string mcmc_message = "";
  NumericVector mcmc_values(3, NA_REAL);
  VectorFactors *factors = new VectorFactors(likelihood_estimators);
  if (mcmc_termination_method.size() > 0) {
    try {
      MCMC *the_mcmc = make_mcmc(model, parameters, &the_data, mcmc_termination_method,
                                 mcmc_weights_method, full_index);
      // the kernels that need transforms or gradients of the particle (as
      // the algorithms pass them)
      std::vector<const ProposalKernel *> proposals = the_mcmc->get_proposals();
      Particle start(point.deep_copy(), factors, &proposals, &proposals);
      start.evaluate_likelihoods(full_index);
      auto one_move = [&](MCMC *mcmc) {
        Particle moved = start;
        Particle scratch;
        rng.seed(seed_value);
        seed = seed_value;
        mcmc->move_in_place(rng, moved, scratch);
        return moved.get_cached_likelihoods(full_index) + arma::accu(moved.parameters.get_colvec(moved.parameters.get_vector_variables()));
      };
      mcmc_values[0] = one_move(the_mcmc);
      MCMC *duplicate = the_mcmc->mcmc_duplicate();
      mcmc_values[1] = one_move(duplicate);
      delete duplicate;
      mcmc_values[2] = one_move(the_mcmc);
      delete the_mcmc;
    } catch (std::exception &e) {
      mcmc_message = e.what();
    }
  }
  delete factors;
  if (full_index != NULL)
    delete full_index;

  return Rcpp::List::create(Named("original") = original, Named("copy") = copy,
                      Named("original_again") = original_again,
                      Named("messages") = messages,
                      Named("mcmc") = mcmc_values,
                      Named("mcmc_message") = mcmc_message);
}

// Record R's main thread when the package is loaded, for the interrupt
// checks (check_interrupt in chunk_parallel.h; dev/plan.md 3.11).
// [[Rcpp::init]]
void ilike_record_main_thread(DllInfo* dll) {
  ilike::ilike_main_thread() = std::this_thread::get_id();
}
