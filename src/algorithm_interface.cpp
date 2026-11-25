#include <algorithm>
#include "vector_factor_variables.h"
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


// (dev/plan.md 3.20, S3) an algorithm whose particles vary a parameter of an
// unnormalised distribution's normaliser would weight them without Z
static void refuse_unknown_normaliser(const RecipeValue &model, const std::string &algorithm)
{
  if (!model.contains("unnormalised_parameters"))
    return;
  std::vector<std::string> names = model["unnormalised_parameters"];
  std::string list;
  for (size_t i = 0; i < names.size(); ++i)
    list += (i ? ", " : "") + names[i];
  stop(algorithm + ": the parameter(s) " + list + " of an unnormalised distribution's normalising constant vary between "
       "particles, and the constant is not known, so the weights would be wrong. Fix them (as data or constants), "
       "or use MCMC with moves that do not change them; doubly intractable methods are planned (dev/plan.md 3.20).");
}

ImportanceSampler *get_importance_sampler(
    RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model,
    const RecipeValue &parameters, const Parameters &algorithm_parameters,
    size_t number_of_importance_points, bool parallel_in, size_t grain_size_in,
    const std::string &results_name_in, size_t *seed,
    std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators,
    const NestedSequence *nested) {
  // (nested: an inner estimator of a factor, integrate u with is(…),
  // dev/plan.md 3.21 step (b): no message, no stored points, and the outer
  // sequences that change the factor in place of the inner model's own)
  refuse_unknown_normaliser(model, "Importance sampling");
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
  if (nested) {
    sequencer_types = nested->types;
    sequencer_variables = nested->variables;
    sequencer_schedules.assign(sequencer_types.size(), std::vector<double>());
  }

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() >= 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  CancellationPlan plan;

  if (model.contains("is_proposal")) {

    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
  } else {
    if (!nested)
      inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, prior_as_proposal(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
  }


  ImportanceSampler *sampler = new ImportanceSampler(
      rng, seed, the_data, algorithm_parameters, number_of_importance_points,
      "", likelihood_estimators, proposal_in, without_cancelled_index, proposal_is_evaluated_in, !nested,
      true, true, false, parallel_in, grain_size_in, "");
  if (nested && !nested->types.empty())
    sampler->set_follows_caller_target(true);
  return sampler;
}

// Output options set by the R wrappers as elements of the recipe.
static void apply_output_options(const RecipeValue &model, SMC *alg, const std::string &results_name) {
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
  // (dev/plan.md 3.21 (c5)) nothing is written without a results name: the
  // SMCs (this one and those inside its factors) then store only the final
  // states of their moves' chains
  alg->set_outputs_written(!results_name.empty());
  // (dev/plan.md 3.21 (c2')) SMC_with_MCMC(..., move_only_after_resampling = TRUE)
  if (model.contains("move_only_after_resampling") && (model["move_only_after_resampling"]).as_bool())
  {
    SMCMCMCMove* moves = dynamic_cast<SMCMCMCMove*>(alg);
    if (moves != NULL)
      moves->set_move_only_after_resampling(true);
  }
}

// After a run with parallel = TRUE (dev/plan.md 2.6a, G4).
static void warn_if_parallel_refused(const SMC *alg) {
  if (alg->parallel_was_refused())
    Rcpp::warning("parallel = TRUE: part of this run was done sequentially, because "
                  "it was not on the chunked path or a factor estimates its likelihood "
                  "with an inner algorithm. The results are the same as a sequential run.");
}

// After an MCMC run: divergent nuts transitions after warm-up (3.18 item 9,
// C212), as Stan warns. A divergence means the leapfrog integrator broke down
// (the Hamiltonian error passed 1000), typically where the target curves
// sharply: the chain may not explore that region, so the estimates may be
// biased even though each transition is valid.
static void warn_if_divergent(const SMCOutput *output) {
  size_t divergent = 0, iterations = 0;
  output->count_divergences(divergent, iterations);
  if (divergent > 0)
    Rcpp::warning("hmc or nuts: " + std::to_string(divergent) + " of " + std::to_string(iterations) +
                  " transitions after warm-up were divergent. The chain may miss regions of high "
                  "curvature: try a smaller step, a transform, or reparameterising the model. "
                  "See load_mcmc_diagnostics().");
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

  apply_output_options(model, alg, results_name_in);
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

  apply_output_options(model, alg, results_name_in);

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
  if (sequencer_info.length() >= 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  CancellationPlan plan;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, plan,
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

  apply_output_options(model, alg, results_name_in);
  alg->set_chains_only(true); // (MCMC: chains, not particles to tune from; C224)

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run();
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);
  warn_if_divergent(output);

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
  if (sequencer_info.length() >= 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  CancellationPlan plan;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, plan,
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

  apply_output_options(model, alg, results_name_in);
  alg->set_chains_only(true); // (MCMC: chains, not particles to tune from; C224)

  // owned here, so that an interrupt (or an error) frees them (dev/plan.md 3.11)
  std::unique_ptr<SMC> alg_owner(alg);
  SMCOutput *output = alg->run(list_to_parameters(fixed_parameter_list));
  std::unique_ptr<SMCOutput> output_owner(output);

  warn_if_parallel_refused(alg);
  warn_if_divergent(output);

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
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators,
    const NestedSequence *nested) {
  // (nested: the inner SMC of integrate u with smc(…), dev/plan.md 3.21 (c1):
  // no message; with outer sequences, (c2): it follows the outer target, its
  // sequence variable the outer one, its factors changing with it)
  refuse_unknown_normaliser(model, "SMC");
  std::vector<LikelihoodEstimator *> likelihood_estimators;
  IndependentProposalKernel *proposal_in;
  bool proposal_is_evaluated_in;
  const bool follows = (nested != NULL) && !nested->types.empty();

  RecipeValue sequencer_info;
  if (follows) {
    if ((nested->types.size() != 1) || !iequals(nested->types[0], "parameter"))
      stop("integrate ... with smc(...): the outer SMC's sequence must be one parameter sequence (e.g. "
           "parameter(epsilon, ...)) for the inner SMC to follow it.");
    // the outer variable, with a placeholder schedule (each output's targets
    // are set as it follows the outer ones: SMCMCMCMove::follow_from_start)
    sequencer_info = RecipeValue::list({std::vector<std::string>{"parameter"},
                                        std::vector<std::string>{nested->variables[0]},
                                        std::vector<std::vector<double>>{std::vector<double>{arma::datum::inf}},
                                        size_t(100)});
  } else {
    sequencer_info = get_smc_sequencer_info(model, parameters, smc_sequencer_method);
  }
  // (dev/plan.md 3.20, S2) data tempering's options
  DataTemperingOptions data_tempering;
  if (sequencer_info.length() >= 6) {
    std::vector<double> options = sequencer_info[5];
    data_tempering.random_order = (options[0] != 0.0);
    data_tempering.anneal = (options[1] != 0.0);
  }
  std::vector<std::string> sequencer_types = sequencer_info[0];
  std::vector<std::string> sequencer_variables = sequencer_info[1];
  std::vector<std::vector<double>> sequencer_schedules = sequencer_info[2];

  std::vector<int> factors_affected_by_smc_sequence;
  if (sequencer_info.length() >= 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;

  CancellationPlan plan;

  if (model.contains("is_proposal")) {
    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators, data_tempering);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
  } else {
    if (nested == NULL)
      inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, prior_as_proposal(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        without_cancelled_index, without_priors_index, full_index,
        plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators, data_tempering);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
  }

  // note that this will always result in the temperature being adapted in the
  // SMC algorithm (not smcfixed)


  MCMC *the_mcmc =
      make_mcmc(model, parameters, the_data, mcmc_termination_method,
                mcmc_weights_method, full_index);

  SMCCriterion *resampling_criterion =
      get_resampling_method(model, parameters, number_of_particles);

  // (dev/plan.md 3.20, S2) with whole data points, the criterion is a step
  // function of the fraction, flat between data points: a bisection converges
  // onto a step from the side where no data point is added, and the sampler
  // makes ever smaller steps without progress. Annealing each data point in
  // makes it continuous.
  if ((sequencer_info.length() >= 6) && !data_tempering.anneal &&
      adaptive_target_method.contains("method") &&
      (adaptive_target_method["method"].as_string() != "positive"))
    stop("data(..., anneal = false) adds whole data points, between which an adaptive target (target = ...) "
         "cannot choose; use anneal = true (each data point annealed in), or a fixed schedule without target.");

  // (dev/plan.md 3.21 (c2''), F1) SMC^2 with each observation annealed in: a
  // step may take in a share of the next observation but not pass a whole
  // one (that would need the filter's propagation, which draws random
  // numbers, inside a reweighting). With an adaptive target each whole
  // observation j/n becomes a point of the schedule, so the search stops
  // there (the sequencer adapts between consecutive points); a fixed
  // schedule must contain them
  if (data_tempering.anneal && !plan.followed_observations.empty() && !sequencer_schedules.empty()) {
    std::vector<double> &schedule = sequencer_schedules.back();
    const bool adaptive = adaptive_target_method.contains("method") &&
                          (adaptive_target_method["method"].as_string() != "positive");
    if (schedule.size() >= 1) {
      const double lo = *std::min_element(schedule.begin(), schedule.end());
      const double hi = *std::max_element(schedule.begin(), schedule.end());
      std::vector<double> points;
      for (size_t n : plan.followed_observations)
        for (size_t j = 0; j <= n; ++j) {
          double point = double(j)/double(n);
          if ((point > lo + 1e-9) && (point < hi - 1e-9))
            points.push_back(point);
        }
      auto in_schedule = [&](double point) {
        for (double value : schedule)
          if (std::abs(value - point) <= 1e-9)
            return true;
        return false;
      };
      if (adaptive) {
        const bool increasing = (schedule.front() <= schedule.back());
        for (double point : points)
          if (!in_schedule(point))
            schedule.push_back(point);
        std::sort(schedule.begin(), schedule.end());
        if (!increasing)
          std::reverse(schedule.begin(), schedule.end());
      } else {
        for (double point : points)
          if (!in_schedule(point))
            stop("SMC^2: data(..., anneal = true) with a particle filter takes in at most one observation per step: a fixed "
                 "schedule must contain every whole observation j/n (missing: " + std::to_string(point) +
                 "); or give an adaptive target (target = ...), which stops at each of them.");
      }
    }
  }

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
  if (follows)
    alg->set_follows_caller_target(nested->variables[0]);
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

  apply_output_options(model, alg, results_name_in);

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

  apply_output_options(model, alg, results_name_in);

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
  if (sequencer_info.length() >= 5) {
    std::vector<int> factors_affected_by_smc_sequence_temp = sequencer_info[4];
    for (auto &k : factors_affected_by_smc_sequence_temp) k = k - 1;
    factors_affected_by_smc_sequence = factors_affected_by_smc_sequence_temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;

  CancellationPlan plan;

  proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

  // the scale the ensemble moves on: transform = … in the algorithm section
  // (C209; the ensemble starts from the transformed prior sample, the
  // simulators get the model's scale, outputs are on the model's scale)
  std::shared_ptr<Transform> transform = NULL;
  if (model.contains("algorithm_transform"))
    transform.reset(transform_from_recipe(model["algorithm_transform"]));

  std::vector<size_t> enk_indices =
      get_enk_likelihood_indices(model, parameters);

  likelihood_estimators = get_likelihood_estimators(
      rng, seed, the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
      without_cancelled_index, without_priors_index, full_index, plan,
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

// MFDS's transform maps the ensemble's (Gaussian) scale to the model's; the
// other direction is not needed (the initial ensemble is given on the
// ensemble's scale)
static Parameters mfds_no_transform_to_ensemble_scale(const Parameters &parameters)
{
  Rcpp::stop("MFDS: internal: no map from the model's scale to the ensemble's.");
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

  // what the lowering of a file in the ilike language prepared for MFDS
  // (dev/plan.md 3.16): a joint simulator of the likelihood, the Gaussian
  // prior's moments, the map to the model's scale, or why they are missing
  RecipeValue mfds;
  bool has_mfds = model.contains("mfds");
  if (has_mfds)
    mfds = model["mfds"];

  // the Gaussian prior (on the scale the ensemble moves on): its means are
  // the artificial observation of the parameters
  Parameters prior_means;
  Parameters prior_covariances;
  std::vector<std::string> log_variables;
  if (has_mfds && mfds.contains("prior")) {
    Data moments = mfds["prior"].function_pointer<DataPtr>()();
    std::vector<std::string> names = moments.get_vector_variables();
    for (const std::string &name : names) {
      if (name.find(':') != std::string::npos)
        continue;
      prior_means[name] = moments[name];
      prior_covariances[name] = moments[name + ":cov"];
      if (std::find(names.begin(), names.end(), name + ":log") != names.end())
        log_variables.push_back(name);
    }
  } else {
    if (has_mfds && mfds.contains("prior_reason"))
      stop("MFDS needs a Gaussian prior (or a log-normal one, or a transformation of a Gaussian), but " +
           mfds["prior_reason"].as_string() + ".");
    get_gaussian_prior_moments(model, parameters, prior_means, prior_covariances, log_variables);
  }

  // the simulator of the data: generated for the whole likelihood, or the
  // one factor with a simulator
  SimulateModelPtr simulate_data_model = NULL;
  std::vector<std::string> data_variables;
  std::shared_ptr<Transform> summary_statistics = NULL;
  if (has_mfds && mfds.contains("simulate")) {
    simulate_data_model = mfds["simulate"].function_pointer<SimulateModelPtr>();
    data_variables = split(mfds["data_variable"].function_pointer<StringPtr>()(), ';');
  } else {
    if (has_mfds && mfds.contains("simulate_reason"))
      stop("MFDS needs a likelihood it can simulate, but " + mfds["simulate_reason"].as_string() + ".");
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
  }
  if (simulate_data_model == NULL)
    stop("MFDS: the model needs a simulate_data_model (and data_variable) for the likelihood.");

  // the map from the Gaussian scale to the model's parameters: exp for
  // log-normal priors, or the user's transform (with which the prior
  // statements give the Gaussian itself)
  TransformPtr to_model = NULL;
  if (model.contains("algorithm_transform"))
    stop("MFDS: transform = auto (or a list) is not available: the prior must be Gaussian on the scale the ensemble "
         "moves on, which log-normal priors already are on the log scale; for other cases give your own transformation "
         "(transform = cpp::f, from the Gaussian scale to the model's).");
  if (model.contains("mfds_transform")) {
    if (!log_variables.empty())
      stop("MFDS: with transform = …, the prior gives the Gaussian on the scale the ensemble moves on, so it must be "
           "norm or mvnorm (" + log_variables[0] + " has a log-normal prior).");
    to_model = model["mfds_transform"].function_pointer<TransformPtr>();
  } else if (!log_variables.empty()) {
    if (!has_mfds || !mfds.contains("log_transform"))
      stop("MFDS: log-normal priors need a recipe in the ilike language.");
    to_model = mfds["log_transform"].function_pointer<TransformPtr>();
  }
  std::shared_ptr<Transform> transform = NULL;
  if (to_model != NULL)
    transform = std::make_shared<Transform>(&mfds_no_transform_to_ensemble_scale, to_model);

  // the initial ensemble, drawn from the prior (with the streams' convention
  // of the ensemble workers: seed, then advance the seed by the ensemble size),
  // on the Gaussian scale (the log of log-normal variables)
  std::unique_ptr<IndependentProposalKernel> prior(get_prior_as_simulate_only_proposal(model, parameters));
  rng.seed(seed);
  std::vector<Parameters> initial_points;
  initial_points.reserve(number_of_ensemble_members);
  for (size_t i = 0; i < number_of_ensemble_members; ++i) {
    Parameters point = prior->independent_simulate(rng);
    for (const std::string &name : log_variables)
      point[name] = arma::log(point[name]);
    initial_points.push_back(point);
  }
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
      summary_statistics, transform, data_variables, false, 100000, results_name_in));

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

  CancellationPlan plan;

  if (model.contains("is_proposal")) {

    proposal_in = get_proposal(model, parameters, the_data);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        full_index, plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
  } else {

    inform_once("No importance sampling proposal specified; using prior.");

    proposal_in = get_prior_as_simulate_only_proposal(model, parameters);

    likelihood_estimators = get_likelihood_estimators(
        rng, seed, the_data, model, parameters, algorithm_parameters, prior_as_proposal(),
        sequencer_types, sequencer_variables, sequencer_schedules, proposal_in,
        evaluated_in_initial_weight_update, evaluated_in_pf_weight_update,
        full_index, plan, factors_affected_by_smc_sequence,
        data_created_in_get_likelihood_estimators,
        data_created_in_get_measurement_covariance_estimators);

    proposal_is_evaluated_in = plan.proposal_is_evaluated; // (dev/plan.md 3.20)
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

  // the bootstrap filter: a transition proposal that is the transition model
  // simulates exactly from it, so both cancel in the weights (dev/plan.md
  // 3.20, rule ii)
  bool transition_proposal_is_evaluated_in;
  if (transition_proposal == NULL) {
    inform_once("No transition proposal specified; using transition model.");

    transition_proposal_is_evaluated_in = kernel_density_is_evaluated(transition_model_as_proposal());

    transition_proposal = transition_model->proposal_kernel_duplicate();
  } else {
    transition_proposal_is_evaluated_in = kernel_density_is_evaluated(KernelDeclaration());

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

  apply_output_options(model, alg, results_name_in);

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

  apply_output_options(model, alg, results_name_in);

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

// The gradient of the log target (the sum of the factors' gradients) with
// respect to one variable at a point: for checking gradients, e.g. those
// generated from model statements (dev/plan.md I4e), against finite differences.
// [[Rcpp::export]]
arma::mat do_evaluate_gradient_log(const Rcpp::List &model_R,
                                   const Rcpp::List &model_parameters_R,
                                   const Rcpp::List &parameters_in,
                                   const std::string &variable,
                                   const IntegerVector &index_in)
{
  RecipeValue model = recipe_from_R(model_R);
  RecipeValue model_parameters = recipe_from_R(model_parameters_R);
  Parameters algorithm_parameters;
  Parameters parameters = list_to_parameters(parameters_in);
  RandomNumberGenerator rng;
  size_t seed = 1;
  Data the_data = get_data(model);
  std::vector<Data> data_created_in_get_likelihood_estimators;
  std::vector<Data> data_created_in_get_measurement_covariance_estimators;
  std::vector<std::string> sequencer_types, sequencer_variables;
  std::vector<std::vector<double>> sequencer_schedules;
  std::vector<int> factors_affected_by_smc_sequence;
  VectorIndex *without_cancelled_index = NULL, *without_priors_index = NULL, *full_index = NULL;
  CancellationPlan plan;
  std::vector<LikelihoodEstimator *> likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, model_parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, plan,
      factors_affected_by_smc_sequence,
      data_created_in_get_likelihood_estimators,
      data_created_in_get_measurement_covariance_estimators);
  std::unique_ptr<Factors> factors(new VectorFactors(likelihood_estimators));
  std::unique_ptr<Index> index((index_in.size() == 0) ? full_index->duplicate()
                                                      : user_factor_index(index_in, likelihood_estimators.size()));
  factors->setup(parameters);
  Particle particle(parameters, factors.get(), NULL, NULL);
  particle.evaluate_smcfixed_part_of_likelihoods(index.get());
  particle.evaluate_smcadaptive_part_given_smcfixed_likelihoods(index.get());
  VectorFactorVariables *factor_variables = dynamic_cast<VectorFactorVariables*>(particle.factor_variables.get());
  if (factor_variables == NULL)
    Rcpp::stop("evaluate_gradient_log: not available for this model.");
  arma::mat gradient;
  if (is_stacked_name(variable))
  {
    // several variables ("a;b;c"): the stacked gradient, as proposals ask for
    // it (each factor asked once for all of them, C204)
    gradient = particle.direct_get_gradient_of_log(variable, index.get());
  }
  else
  {
    gradient.zeros(arma::size(static_cast<const Parameters&>(parameters)[variable]));
    for (auto i = index->begin(); i != index->end(); ++i)
    {
      LikelihoodEstimatorOutput *output = factor_variables->get_output(*i);
      if (output != NULL)
        gradient = gradient + output->get_gradient_of_log(variable, particle.parameters);
    }
  }
  delete without_cancelled_index;
  delete without_priors_index;
  delete full_index;
  return gradient;
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
  CancellationPlan plan;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, model_parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, plan,
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
  CancellationPlan plan;

  std::vector<LikelihoodEstimator *> likelihood_estimators;
  likelihood_estimators = get_likelihood_estimators(
                                                    &rng, &seed, &the_data, model, model_parameters, algorithm_parameters, KernelDeclaration(),
                                                    sequencer_types, sequencer_variables, sequencer_schedules, NULL,
                                                    without_cancelled_index, without_priors_index, full_index, plan,
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
  if (sequencer_info.length() >= 5) {
    std::vector<int> temp = sequencer_info[4];
    for (auto &k : temp) k = k - 1;
    factors_affected_by_smc_sequence = temp;
  }

  VectorIndex *without_cancelled_index = NULL;
  VectorIndex *without_priors_index = NULL;
  VectorIndex *full_index = NULL;
  CancellationPlan plan;
  std::vector<LikelihoodEstimator *> likelihood_estimators = get_likelihood_estimators(
      &rng, &seed, &the_data, model, parameters, algorithm_parameters, KernelDeclaration(),
      sequencer_types, sequencer_variables, sequencer_schedules, NULL,
      without_cancelled_index, without_priors_index, full_index, plan,
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


// Draws from the simulator of an unnormalised distribution (dev/plan.md
// 3.20, S3): method "exact" (coupling from the past) or "sweeps" (that many
// heat-bath sweeps from the value of the variable in parameters, or from the
// base measure if it is not there). Internal: for tests and later methods.
typedef void (*UnnormalisedSweepPtr)(RandomNumberGenerator &rng, Parameters &parameters, const Data &data);
typedef Parameters (*UnnormalisedExactPtr)(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data);

// [[Rcpp::export]]
arma::mat do_simulate_unnormalised(const Rcpp::List &model_R,
                                   const std::string &variable,
                                   const Rcpp::List &parameters_in,
                                   size_t draws,
                                   const std::string &method,
                                   size_t sweeps,
                                   size_t seed)
{
  RecipeValue model = recipe_from_R(model_R);
  if (!model.contains("unnormalised"))
    stop("The recipe has no unnormalised distribution.");
  const RecipeValue &entries = model["unnormalised"];
  int found = -1;
  for (size_t i = 0; i < entries.size(); ++i)
    if (entries[i]["variable"].as_string() == variable)
      found = int(i);
  if (found < 0)
    stop("The recipe has no unnormalised distribution for " + variable + ".");
  const RecipeValue &entry = entries[size_t(found)];
  Data the_data = get_data(model);
  Parameters start = list_to_parameters(parameters_in);
  RandomNumberGenerator rng;
  rng.seed(seed);
  arma::mat output;
  for (size_t d = 0; d < draws; ++d)
  {
    arma::mat x;
    if (method == "exact")
    {
      if (!entry.contains("exact"))
        stop("No exact simulator for " + variable + " (" + entry["kind"].as_string() + ").");
      Parameters draw = entry["exact"].function_pointer<UnnormalisedExactPtr>()(rng, start, the_data);
      x = draw[variable];
    }
    else
    {
      Parameters state = start.deep_copy();
      UnnormalisedSweepPtr sweep = entry["sweep"].function_pointer<UnnormalisedSweepPtr>();
      for (size_t k = 0; k < std::max(sweeps, size_t(1)); ++k)
        sweep(rng, state, the_data);
      x = state[variable];
    }
    if (output.n_rows == 0)
      output.set_size(draws, x.n_elem);
    output.row(d) = arma::vectorise(x).t();
  }
  return output;
}
