#include "recipe_internal.h"
#include <set>
#include "plate_likelihood_estimator.h"

// The target of an algorithm, from the recipe's factors and its SMC sequence
// (dev/plan.md 3.20, S1): each factor's kind, whether the sequence names it,
// and whether the sequence is annealing. plan_cancellations
// (cancellation_plan.cpp) then decides what cancels.
TargetSpecification target_specification(
    const RecipeValue &model, const std::vector<std::string> &sequencer_types,
    const std::vector<std::string> &sequencer_variables,
    const std::vector<int> &factors_affected_by_smc_sequence) {
  TargetSpecification target;
  for (size_t i = 0; i < sequencer_types.size(); ++i) {
    if ((iequals(sequencer_types[i], "annealing")) ||
        (iequals(sequencer_types[i], "tempering"))) {
      // (the first annealing sequence)
      target.annealing = true;
      target.annealing_variable = sequencer_variables[i];
      break;
    }
  }
  target.sequence = !sequencer_types.empty();
  target.number_of_sequences = sequencer_types.size();
  for (size_t i = 0; i < sequencer_types.size(); ++i)
    if (iequals(sequencer_types[i], "data")) {
      target.data_tempering = true;
      target.data_variable = sequencer_variables[i];
    }
  target.factors_listed = !factors_affected_by_smc_sequence.empty();
  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];
    for (size_t j = 0; j < factors_affected_by_smc_sequence.size(); ++j) {
      if ((factors_affected_by_smc_sequence[j] < 0) ||
          (size_t(factors_affected_by_smc_sequence[j]) >= factors.size()))
        stop("smc_sequence: the factors affected must be between 1 and " +
             std::to_string(factors.size()) + ".");
    }
    for (size_t i = 0; i < factors.size(); ++i) {
      bool listed = false;
      for (size_t j = 0; j < factors_affected_by_smc_sequence.size(); ++j)
        if (int(i) == factors_affected_by_smc_sequence[j])
          listed = true;
      target.listed.push_back(listed);
      // (in the order get_likelihood_estimators reads the entry; a malformed
      // entry is reported there)
      FactorKind kind = FactorKind::likelihood;
      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];
        if (current_factor.contains("prior") ||
            current_factor.contains("evaluate_log_prior"))
          kind = FactorKind::prior;
        else if (current_factor.contains("simulate_prior"))
          kind = FactorKind::none;
      }
      target.kinds.push_back(kind);
    }
  }
  return target;
}

// A factor written as model statements (I4b): a plate of terms, when its
// entry has "plate"; otherwise the exact estimator of its whole function.
LikelihoodEstimator *make_model_statement_estimator(
    const RecipeValue &factor, const RecipeValue &model,
    RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in,
    bool smcfixed) {
  const RecipeValue &plate = factor["plate"];
  std::vector<std::string> variables;
  const RecipeValue &names = plate["variables"];
  for (size_t k = 0; k < names.size(); ++k)
    variables.push_back(names[k].as_string());
  bool check = model.contains("check_factor_index") &&
               model["check_factor_index"].as_bool();
  PlateLikelihoodEstimator *estimator = new PlateLikelihoodEstimator(
      rng_in, seed_in, data_in,
      plate["structure"].function_pointer<PlateStructurePtr>(),
      plate["term"].function_pointer<PlateTermPtr>(), variables,
      plate["name"].as_string(), check, smcfixed);
  // its gradient, generated from its terms (I4e)
  if (plate.contains("gradient"))
    estimator->set_gradient(plate["gradient"].function_pointer<PlateGradientPtr>());
  if (plate.contains("gradients"))
    estimator->set_gradients(plate["gradients"].function_pointer<PlateGradientsPtr>());
  if (plate.contains("value_gradients"))
    estimator->set_value_gradients(plate["value_gradients"].function_pointer<PlateValueGradientsPtr>());
  return estimator;
}

std::vector<LikelihoodEstimator *> get_likelihood_estimators(
    RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in,
    const RecipeValue &model, const RecipeValue &model_parameters,
    const Parameters &algorithm_parameters,
    const KernelDeclaration &importance_proposal,
    const std::vector<std::string> &sequencer_types,
    const std::vector<std::string> &sequencer_variables,
    const std::vector<std::vector<double>> &sequencer_schedules,
    IndependentProposalKernel *proposal_in,
    VectorIndex *&without_cancelled_index, VectorIndex *&without_priors_index,
    VectorIndex *&full_index, CancellationPlan &plan,
    const std::vector<int> &factors_affected_by_smc_sequence,
    std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators,
    const DataTemperingOptions &data_tempering) {

  std::vector<size_t> without_cancelled_index_vector;
  std::vector<size_t> without_priors_index_vector;
  std::vector<size_t> full_index_vector;

  // what cancels (dev/plan.md 3.20, S1)
  plan = plan_cancellations(
      target_specification(model, sequencer_types, sequencer_variables,
                           factors_affected_by_smc_sequence),
      importance_proposal);
  const std::string &annealing_variable = plan.annealing_variable;

  PowerFunctionPtr power = annealing_power;
  PowerFunctionPtr second_power = annealing_one_minus_power;

  std::vector<LikelihoodEstimator *> likelihood_estimators;

  if (model.contains("factor")) {

    RecipeValue factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {

      // (the plan: annealed, fixed along the sequence, in the weights)
      const FactorPlan &factor_plan = plan.factors[i];
      const bool factor_is_annealed = factor_plan.annealed;
      const bool factor_is_fixed_in_smc = factor_plan.fixed_in_smc;

      // (S2) data tempering takes the terms of a factor written as model
      // statements one at a time; another factor cannot be tempered this way
      // (dev/plan.md 3.21 (c2')) a nested particle filter takes its
      // observations in one at a time (SMC^2)
      // (3.21 (c2''')) and a Kalman or ensemble Kalman filter
      const bool tempered_filter = factor_plan.data_tempered && (factors[i]).is_list() &&
        factors[i].contains("algorithmic_likelihood") && factors[i]["algorithmic_likelihood"].contains("type") &&
        (iequals(factors[i]["algorithmic_likelihood"]["type"].as_string(), "pf") ||
         iequals(factors[i]["algorithmic_likelihood"]["type"].as_string(), "kf") ||
         iequals(factors[i]["algorithmic_likelihood"]["type"].as_string(), "enkf"));
      if (factor_plan.data_tempered && !tempered_filter &&
          !((factors[i]).is_list() && factors[i].contains("evaluate_log_likelihood") &&
            factors[i].contains("plate")))
        stop("smc_sequence data tempering: factor " + std::to_string(i + 1) +
             " is not written as model statements, so its data points cannot be taken in one at a time. "
             "List the factors to temper (data(..., factors = [...])), or read the sequence variable in the "
             "factor's own function with parameter(...).");

      if ((factors[i]).is_list()) {

        RecipeValue current_factor = factors[i];
        if (current_factor.contains("prior")) {
          DistributionFactor *new_factor;

          if ((current_factor["prior"]).is_list()) {

            RecipeValue current_prior = current_factor["prior"];
            if (current_prior.contains("type") &&
                current_prior.contains("variables")) {
              new_factor = make_builtin_prior(current_prior, model_parameters);

              LikelihoodEstimator *new_likelihood_estimator =
                  new ExactLikelihoodEstimator(rng_in, seed_in, data_in,
                                               new_factor, true);

              full_index_vector.push_back(likelihood_estimators.size());
              if (factor_plan.in_weight)
                without_cancelled_index_vector.push_back(
                    likelihood_estimators.size());
              if (factor_is_annealed) {
                likelihood_estimators.push_back(
                    new AnnealedLikelihoodEstimator(
                        rng_in, seed_in, data_in, new_likelihood_estimator,
                        power, annealing_variable, false));
              } else {
                likelihood_estimators.push_back(new_likelihood_estimator);
              }
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.contains("evaluate_log_prior")) {
          RecipeValue evaluate_log_prior_value = current_factor["evaluate_log_prior"];
          // The gradient is optional: it is only needed by gradient-based
          // proposals (e.g. Langevin, Barker, HMC).
          DistributionFactor *new_factor;
          if (current_factor.contains(
                  "evaluate_gradient_log_prior")) {
            RecipeValue evaluate_gradient_log_prior_value =
                current_factor["evaluate_gradient_log_prior"];
            new_factor = new CustomDistributionFactor(
                (evaluate_log_prior_value).function_pointer<EvaluateLogDistributionPtr>(),
                (evaluate_gradient_log_prior_value).function_pointer<EvaluateGradientLogDistributionPtr>());
          } else {
            new_factor = new CustomDistributionFactor(
                (evaluate_log_prior_value).function_pointer<EvaluateLogDistributionPtr>());
          }
          LikelihoodEstimator *new_likelihood_estimator;
          if (current_factor.contains("plate")) {
            delete new_factor;
            new_likelihood_estimator = make_model_statement_estimator(
                current_factor, model, rng_in, seed_in, data_in,
                factor_is_fixed_in_smc);
          } else {
            new_likelihood_estimator =
                new ExactLikelihoodEstimator(rng_in, seed_in, data_in, new_factor,
                                             factor_is_fixed_in_smc);
          }
          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.contains("simulate_prior")) {
        } else if (current_factor.contains(
                       "evaluate_log_likelihood")) {

          RecipeValue evaluate_log_likelihood_value =
              current_factor["evaluate_log_likelihood"];
          // The gradient is optional: it is only needed by gradient-based
          // proposals (e.g. Langevin, Barker, HMC).
          LikelihoodFactor *new_factor;
          if (current_factor.contains(
                  "evaluate_gradient_log_likelihood")) {
            RecipeValue evaluate_gradient_log_likelihood_value =
                current_factor["evaluate_gradient_log_likelihood"];
            new_factor = new CustomLikelihoodFactor(
                (evaluate_log_likelihood_value).function_pointer<EvaluateLogLikelihoodPtr>(),
                (evaluate_gradient_log_likelihood_value).function_pointer<EvaluateGradientLogLikelihoodPtr>(),
                data_in);
          } else {
            new_factor = new CustomLikelihoodFactor(
                (evaluate_log_likelihood_value).function_pointer<EvaluateLogLikelihoodPtr>(),
                data_in);
          }
          // (dev/plan.md 3.21) density(…, maximum = cpp::f): its largest value
          // over the state, for early stopping in particle filters
          if (current_factor.contains("maximum_log_likelihood"))
            static_cast<CustomLikelihoodFactor*>(new_factor)->set_maximum(
                current_factor["maximum_log_likelihood"].function_pointer<EvaluateLogLikelihoodPtr>());
          LikelihoodEstimator *new_likelihood_estimator;
          if (current_factor.contains("plate")) {
            delete new_factor;
            new_likelihood_estimator = make_model_statement_estimator(
                current_factor, model, rng_in, seed_in, data_in,
                factor_is_fixed_in_smc);
            // (dev/plan.md 3.21 (c5)) an ABC kernel's summaries
            if (current_factor.contains("abc_summary") && current_factor.contains("abc_from_summary"))
              static_cast<PlateLikelihoodEstimator *>(new_likelihood_estimator)->set_kernel_summary(
                  current_factor["abc_summary"].function_pointer<AbcSummaryPtr>(),
                  current_factor["abc_from_summary"].function_pointer<AbcFromSummaryPtr>());
            // (dev/plan.md 3.21 (c4)) its maximum, generated from its terms
            if (current_factor.contains("maximum_log_likelihood"))
              static_cast<PlateLikelihoodEstimator *>(new_likelihood_estimator)->set_maximum(
                  current_factor["maximum_log_likelihood"].function_pointer<EvaluateLogLikelihoodPtr>());
            // (S2) its data points taken in along the sequence
            if (factor_plan.data_tempered)
              static_cast<PlateLikelihoodEstimator *>(new_likelihood_estimator)
                  ->set_data_tempering(plan.data_variable, data_tempering.anneal,
                                       data_tempering.random_order,
                                       (seed_in != NULL) ? *seed_in : 0);
          } else {
            new_likelihood_estimator =
                new ExactLikelihoodEstimator(rng_in, seed_in, data_in, new_factor,
                                             factor_is_fixed_in_smc);
          }

          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
          if (factor_plan.in_without_priors)
            without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.contains("sbi_likelihood")) {
          SimulateModelPtr simulate_data_model_in;
          if (current_factor.contains("simulate_data_model")) {
            RecipeValue simulate_data_model_value =
                current_factor["simulate_data_model"];
            simulate_data_model_in =
                (simulate_data_model_value).function_pointer<SimulateModelPtr>();
          } else {
            stop(
                "sbi_likelihood factor must also contain simulate_data_model.");
          }

          std::vector<std::string> data_variables_in;
          if (current_factor.contains("data_variable")) {
            RecipeValue data_variable_value = current_factor["data_variable"];
            std::string augmented_data_variable_names =
                (data_variable_value).function_pointer<StringPtr>()();
            data_variables_in = split(augmented_data_variable_names, ';');
          } else {
            stop("sbi_likelihood factor must also contain data_variable.");
          }

          LikelihoodEstimator *new_likelihood_estimator;

          RecipeValue current_sbi = current_factor["sbi_likelihood"];
          if (current_sbi.contains("type") &&
              current_sbi.contains("parameters")) {
            // the registry's name (an error listing the alternatives if unknown)
            std::string type = find_builtin(current_sbi["type"], BuiltinRole::SBILikelihood).name;

            if ((iequals(type, "euclidean_uniform_abc")) ||
                (iequals(type, "lp_uniform_abc"))) {
              size_t number_of_points;
              std::string tolerance_variable;
              double tolerance;
              double p;
              bool store_output;
              bool parallel;
              size_t grain_size;
              std::string scale_variable_in;
              arma::colvec scale_in;

              if (iequals(type, "euclidean_uniform_abc")) {
                RecipeValue info = get_abc_euclidean_uniform_parameter_info(
                    model_parameters, current_sbi, type);

                number_of_points = info[0];
                std::string temp_tolerance_variable = info[1];
                tolerance_variable = temp_tolerance_variable;
                tolerance = info[2];
                p = 2.0;
                store_output = info[3];
                std::string temp_scale_variable_in = info[4];
                scale_variable_in = temp_scale_variable_in;

                arma::colvec temp_scale_in = info[5];
                scale_in = temp_scale_in;
                // a scale given without a variable's name: kept under a name of its own
                // (it used to be ignored: the kernel then used a scale of 1; dev/plan.md C253)
                if (scale_variable_in.empty() && (scale_in.n_elem > 0))
                  scale_variable_in = "ilike_abc_scale";

                parallel = info[6];
                grain_size = info[7];
              } else {
                RecipeValue info = get_abc_lp_uniform_parameter_info(
                    model_parameters, current_sbi, type);

                number_of_points = info[0];
                std::string temp_tolerance_variable = info[1];
                tolerance_variable = temp_tolerance_variable;
                tolerance = info[2];
                p = info[3];
                store_output = info[4];
                std::string temp_scale_variable_in = info[5];
                scale_variable_in = temp_scale_variable_in;
                arma::colvec temp_scale_in = info[6];
                scale_in = temp_scale_in;
                // a scale given without a variable's name: kept under a name of its own
                // (it used to be ignored: the kernel then used a scale of 1; dev/plan.md C253)
                if (scale_variable_in.empty() && (scale_in.n_elem > 0))
                  scale_variable_in = "ilike_abc_scale";

                parallel = info[7];
                grain_size = info[8];
              }

              bool adaptive = false;
              for (auto k = sequencer_variables.begin();
                   k != sequencer_variables.end(); ++k) {
                if (*k == tolerance_variable) {
                  adaptive = true;
                }
              }

              if (current_factor.contains("summary_statistics")) {
                RecipeValue summary_statistics_value =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

                data_created_in_get_likelihood_estimators.push_back(
                    summary_stats(*data_in));

                if (adaptive == false) {
                  new_likelihood_estimator =
                      make_fixed_epsilon_lp_uniform_abc_likelihood(
                          rng_in, seed_in,
                          &data_created_in_get_likelihood_estimators.back(), p,
                          scale_in, data_variables_in, scale_variable_in,
                          tolerance_variable, tolerance, simulate_data_model_in,
                          std::make_shared<Transform>(summary_stats),
                          number_of_points, store_output, parallel, grain_size);
                } else {
                  new_likelihood_estimator =
                      make_varying_epsilon_lp_uniform_abc_likelihood(
                          rng_in, seed_in,
                          &data_created_in_get_likelihood_estimators.back(), p,
                          scale_in, data_variables_in, scale_variable_in,
                          tolerance_variable, simulate_data_model_in,
                          std::make_shared<Transform>(summary_stats),
                          number_of_points, store_output, parallel, grain_size);
                }
              } else {
                if (adaptive == false) {
                  new_likelihood_estimator =
                      make_fixed_epsilon_lp_uniform_abc_likelihood(
                          rng_in, seed_in, data_in, p, scale_in,
                          data_variables_in, scale_variable_in,
                          tolerance_variable, tolerance, simulate_data_model_in,
                          number_of_points, store_output, parallel, grain_size);
                } else {
                  new_likelihood_estimator =
                      make_varying_epsilon_lp_uniform_abc_likelihood(
                          rng_in, seed_in, data_in, p, scale_in,
                          data_variables_in, scale_variable_in,
                          tolerance_variable, simulate_data_model_in,
                          number_of_points, store_output, parallel, grain_size);
                }
              }
            } else if (iequals(type, "gaussian_abc")) {
              size_t number_of_points;
              std::string tolerance_variable;
              double tolerance;
              bool store_output;
              std::string scale_variable_in;
              arma::colvec scale_in;
              bool parallel;
              size_t grain_size;

              // same function as getting a Gaussian
              RecipeValue info = get_abc_euclidean_uniform_parameter_info(
                  model_parameters, current_sbi, type);

              number_of_points = info[0];
              std::string temp_tolerance_variable = info[1];
              tolerance_variable = temp_tolerance_variable;
              tolerance = info[2];
              store_output = info[3];
              std::string temp_scale_variable_in = info[4];
              scale_variable_in = temp_scale_variable_in;
              arma::colvec temp_scale_in = info[5];
              scale_in = temp_scale_in;
              // a scale given without a variable's name: kept under a name of its own
              // (it used to be ignored: the kernel then used a scale of 1; dev/plan.md C253)
              if (scale_variable_in.empty() && (scale_in.n_elem > 0))
                scale_variable_in = "ilike_abc_scale";

              parallel = info[6];
              grain_size = info[7];

              bool adaptive = false;
              for (auto k = sequencer_variables.begin();
                   k != sequencer_variables.end(); ++k) {
                if (*k == tolerance_variable) {
                  adaptive = true;
                }
              }

              if (current_factor.contains("summary_statistics")) {
                RecipeValue summary_statistics_value =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

                data_created_in_get_likelihood_estimators.push_back(
                    summary_stats(*data_in));

                if (adaptive == false) {
                  new_likelihood_estimator =
                      make_fixed_epsilon_gaussian_abc_likelihood(
                          rng_in, seed_in,
                          &data_created_in_get_likelihood_estimators.back(),
                          scale_in, data_variables_in, scale_variable_in,
                          tolerance_variable, tolerance, simulate_data_model_in,
                          std::make_shared<Transform>(summary_stats),
                          number_of_points, store_output, parallel, grain_size);
                } else {
                  new_likelihood_estimator =
                      make_varying_epsilon_gaussian_abc_likelihood(
                          rng_in, seed_in,
                          &data_created_in_get_likelihood_estimators.back(),
                          scale_in, data_variables_in, scale_variable_in,
                          tolerance_variable, simulate_data_model_in,
                          std::make_shared<Transform>(summary_stats),
                          number_of_points, store_output, parallel, grain_size);
                }
              } else {
                if (adaptive == false) {
                  new_likelihood_estimator =
                      make_fixed_epsilon_gaussian_abc_likelihood(
                          rng_in, seed_in, data_in, scale_in, data_variables_in,
                          scale_variable_in, tolerance_variable, tolerance,
                          simulate_data_model_in, number_of_points,
                          store_output, parallel, grain_size);
                } else {
                  new_likelihood_estimator =
                      make_varying_epsilon_gaussian_abc_likelihood(
                          rng_in, seed_in, data_in, scale_in, data_variables_in,
                          scale_variable_in, tolerance_variable,
                          simulate_data_model_in, number_of_points,
                          store_output, parallel, grain_size);
                }
              }
            } else if (iequals(type, "ienki_abc")) {
              size_t number_of_points;
              size_t estimator_type_in;
              std::string tolerance_variable;
              double epsilon;
              size_t enki_lag;
              std::string shifter_name;
              size_t number_of_targets_in;
              bool enki_on_summary;
              double significance_level;
              bool parallel;
              size_t grain_size;

              RecipeValue info = get_abc_enki_parameter_info(model_parameters,
                                                      current_sbi, type);

              number_of_points = info[0];
              estimator_type_in = info[1];
              std::string temp_tolerance_variable = info[2];
              tolerance_variable = temp_tolerance_variable;
              epsilon = info[3];
              std::string temp_shifter_name = info[4];
              shifter_name = temp_shifter_name;
              number_of_targets_in = info[5];
              significance_level = info[6];
              enki_lag = info[7];
              std::string scale_variable_in = info[8];
              arma::colvec scale_in = info[9];
              enki_on_summary = info[10];
              parallel = info[11];
              grain_size = info[12];

              if ((epsilon < 0.0) ||
                  ((epsilon == 0.0) && (estimator_type_in > 2))) {
                stop("Negative epsilon, or zero epsilon with estimator_type "
                     "not equal to 1 or 2, not allowed in IEnKI-ABC.");
              }

              EnsembleShifter *shifter = make_builtin_enk_shifter(shifter_name);

              bool adaptive = false;
              for (auto k = sequencer_variables.begin();
                   k != sequencer_variables.end(); ++k) {
                if (*k == tolerance_variable) {
                  adaptive = true;
                }
              }

              if (current_factor.contains("summary_statistics")) {
                RecipeValue summary_statistics_value =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

                data_created_in_get_likelihood_estimators.push_back(
                    summary_stats(*data_in));

                if (adaptive == false) {
                  // If epsilon is equal to 0, then use the synthetic
                  // likelihood.

                  if ((epsilon == 0.0) && (estimator_type_in == 0)) {
                    new_likelihood_estimator = make_sl_likelihood(
                                                                  rng_in, seed_in,
                                                                  &data_created_in_get_likelihood_estimators.back(),
                                                                  data_variables_in, simulate_data_model_in,
                                                                  std::make_shared<Transform>(summary_stats),
                                                                  number_of_points, false, false, parallel, grain_size);
                  }
                  else if ((epsilon == 0.0) && (estimator_type_in == 1)) {
                    new_likelihood_estimator = make_sl_likelihood(
                        rng_in, seed_in,
                        &data_created_in_get_likelihood_estimators.back(),
                        data_variables_in, simulate_data_model_in,
                        std::make_shared<Transform>(summary_stats),
                        number_of_points, false, false, parallel, grain_size);
                  } else if ((epsilon == 0.0) && (estimator_type_in == 2)) {
                    new_likelihood_estimator = make_sl_likelihood(
                        rng_in, seed_in,
                        &data_created_in_get_likelihood_estimators.back(),
                        data_variables_in, simulate_data_model_in,
                        std::make_shared<Transform>(summary_stats),
                        number_of_points, true, false, parallel, grain_size);
                  } else {
                    new_likelihood_estimator =
                        make_fixed_epsilon_enki_abc_likelihood(
                            rng_in, seed_in,
                            &data_created_in_get_likelihood_estimators.back(),
                            enki_lag, shifter, scale_variable_in,
                            tolerance_variable, epsilon, number_of_targets_in,
                            scale_in, data_variables_in, simulate_data_model_in,
                            std::make_shared<Transform>(summary_stats),
                            enki_on_summary, number_of_points,
                            significance_level, estimator_type_in, parallel,
                            grain_size);
                  }
                } else {
                  stop("Adaptive approach to choosing epsilon not yet "
                       "implemented in EnKI-ABC.");
                }
              } else {
                if (adaptive == false) {
                  if ((epsilon == 0.0) && (estimator_type_in == 0)) {
                    new_likelihood_estimator = make_sl_likelihood(
                                                                  rng_in, seed_in, data_in, data_variables_in,
                                                                  simulate_data_model_in, number_of_points, false, false,
                                                                  parallel, grain_size);
                  }
                  else if ((epsilon == 0.0) && (estimator_type_in == 1)) {
                    new_likelihood_estimator = make_sl_likelihood(
                        rng_in, seed_in, data_in, data_variables_in,
                        simulate_data_model_in, number_of_points, false, false,
                        parallel, grain_size);
                  } else if ((epsilon == 0.0) && (estimator_type_in == 2)) {
                    new_likelihood_estimator = make_sl_likelihood(
                        rng_in, seed_in, data_in, data_variables_in,
                        simulate_data_model_in, number_of_points, true, false,
                        parallel, grain_size);
                  } else {
                    new_likelihood_estimator =
                        make_fixed_epsilon_enki_abc_likelihood(
                            rng_in, seed_in, data_in, enki_lag, shifter,
                            scale_variable_in, tolerance_variable, epsilon,
                            number_of_targets_in, scale_in, data_variables_in,
                            simulate_data_model_in, number_of_points,
                            significance_level, estimator_type_in, parallel,
                            grain_size);
                  }
                } else {
                  stop("Adaptive approach to choosing epsilon not yet "
                       "implemented in EnKI-ABC.");
                }
              }
            } else if (iequals(type, "sl")) {
              size_t number_of_points;
              bool unbiased;
              bool store_output;
              bool parallel;
              size_t grain_size;

              // same function as getting a Gaussian
              RecipeValue info =
                  get_sl_parameter_info(model_parameters, current_sbi, type);

              number_of_points = info[0];
              unbiased = info[1];
              store_output = info[2];
              parallel = info[3];
              grain_size = info[4];

              if (current_factor.contains("summary_statistics")) {
                RecipeValue summary_statistics_value =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

                data_created_in_get_likelihood_estimators.push_back(
                    summary_stats(*data_in));

                new_likelihood_estimator = make_sl_likelihood(
                    rng_in, seed_in,
                    &data_created_in_get_likelihood_estimators.back(),
                    data_variables_in, simulate_data_model_in,
                    std::make_shared<Transform>(summary_stats),
                    number_of_points, unbiased, store_output, parallel,
                    grain_size);
              } else {
                new_likelihood_estimator = make_sl_likelihood(
                    rng_in, seed_in, data_in, data_variables_in,
                    simulate_data_model_in, number_of_points, unbiased,
                    store_output, parallel, grain_size);
              }
            } else {
              Rcout << "SBI type " << type;
              stop("SBI type unknown.");
            }
          } else {
            stop("Missing information for SBI in model file.");
          }

          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
          if (factor_plan.in_without_priors)
            without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (is_linear_gaussian_data_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_data_model(current_factor, model_parameters);
          std::vector<std::string> state_variable_names = lg.state_variables;
          size_t type = lg.is_function ? 2 : 1;
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new LinearGaussianNoiseProposalKernel(
                lg.variable,
                state_variable_names, // lg.state_variable_string,
                lg.matrix,
                lg.covariance);
          } else if (type == 2) {
            proposal = new LinearGaussianNoiseFunctionProposalKernel(
                lg.variable,
                state_variable_names, // lg.state_variable_string,
                lg.matrix_function,
                lg.covariance_function);
          } else {
            stop("linear_gaussian_data specificiation is invalid.");
          }

          LikelihoodEstimator *new_likelihood_estimator =
              new ExactLikelihoodEstimator(rng_in, seed_in, data_in, proposal,
                                           factor_is_fixed_in_smc);

          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
          if (factor_plan.in_without_priors)
            without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.contains(
                       "nonlinear_gaussian_data_covariance")) {
          if (!current_factor.contains(
                  "nonlinear_gaussian_data_function")) {
            stop("nonlinear_gaussian_data must contain both function and "
                 "covariance.");
          }

          if (!current_factor.contains("type")) {
            stop("nonlinear_gaussian_data must contain a type.");
          }

          if (!current_factor.contains(
                  "nonlinear_gaussian_data_variable")) {
            stop("nonlinear_gaussian_data must contain "
                 "nonlinear_gaussian_data_variable.");
          }

          // std::vector<std::string> variables =
          // data_in->get_vector_variables(); if (variables.size()!=1)
          //   stop("For this linear-Gaussian likelihood we can only have one
          //   variable present in the data.");

          RecipeValue nonlinear_gaussian_data_function_value =
              current_factor["nonlinear_gaussian_data_function"];
          RecipeValue nonlinear_gaussian_data_covariance_value =
              current_factor["nonlinear_gaussian_data_covariance"];

          RecipeValue nonlinear_gaussian_data_variable_value =
              current_factor["nonlinear_gaussian_data_variable"];

          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>(
                    (nonlinear_gaussian_data_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_data_covariance_value).function_pointer<GetNoParamsMatrixPtr>()());
          } else if (type == 2) {
            proposal = new NonLinearGaussianNoiseFunctionProposalKernel(
                (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>(
                    (nonlinear_gaussian_data_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_data_covariance_value).function_pointer<GetMatrixPtr>());
          } else {
            stop("nonlinear_gaussian_data specificiation is invalid.");
          }

          LikelihoodEstimator *new_likelihood_estimator =
              new ExactLikelihoodEstimator(rng_in, seed_in, data_in, proposal,
                                           factor_is_fixed_in_smc);

          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
          if (factor_plan.in_without_priors)
            without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.contains(
                       "algorithmic_likelihood")) {
          if ((current_factor["algorithmic_likelihood"]).is_list()) {
            RecipeValue llhd_details = current_factor["algorithmic_likelihood"];

            if (!llhd_details.contains("type")) {
              stop("algorithmic_likelihood must contain a type.");
            }

            if (!llhd_details.contains("model")) {
              stop("algorithmic_likelihood must contain a model.");
            }

            if (!llhd_details.contains("parameters")) {
              stop("algorithmic_likelihood must contain parameters.");
            }

            std::string llhd_type =
                find_builtin(llhd_details["type"], BuiltinRole::AlgorithmicLikelihood).name;
            RecipeValue llhd_model = llhd_details["model"];

            LikelihoodEstimator *new_likelihood_estimator;
            if (iequals(llhd_type, "kf")) {
              RecipeValue kf_params =
                  get_kf_info(model_parameters, llhd_details, llhd_type);

              size_t iterations_to_store = kf_params[0];
              double steady_state = kf_params[1];

              KalmanFilter *inner_kf = get_kalman_filter(
                  data_in, llhd_model, model_parameters,
                  algorithm_parameters, iterations_to_store, false, "");
              // (C276) the steady-state option: not with SMC^2, whose outputs
              // take their steps in turn (the frozen state is the filter's)
              if (steady_state > 0.0)
              {
                if (factor_plan.data_tempered)
                {
                  delete inner_kf;
                  stop("kf: steady_state is not available when the caller's sequence takes in the filter's observations one at a time (SMC^2).");
                }
                inner_kf->set_steady_state_tolerance(steady_state);
              }
              // (dev/plan.md 3.21 (c2''')) SMC^2 with a Kalman filter: the
              // caller's data tempering takes in its observations one at a time
              if (factor_plan.data_tempered)
              {
                inner_kf->set_follows_caller_target(plan.data_variable, data_tempering.anneal);
                plan.followed_observations.push_back(inner_kf->observations());
              }
              new_likelihood_estimator = inner_kf;
            } else if (iequals(llhd_type, "enkf")) {
              RecipeValue enkf_params =
                  get_enkf_info(model_parameters, llhd_details, llhd_type);

              size_t number_of_ensemble_members = enkf_params[0];
              std::string shifter_name = enkf_params[1];
              size_t iterations_to_store = enkf_params[2];
              bool parallel = enkf_params[3];
              size_t grain_size = enkf_params[4];

              // (the square-root shifter needs an "A" (H) matrix built from the
              // state variables in EnsembleKalman's packing instructions and
              // those stored in the measurement covariance estimator, in the
              // setup functions called before running the method)
              EnsembleShifter *shifter = make_builtin_enk_shifter(shifter_name);

              EnsembleKalmanFilter *inner_enkf = get_ensemble_kalman_filter(
                  rng_in, data_in, llhd_model, model_parameters,
                  algorithm_parameters, number_of_ensemble_members, shifter,
                  iterations_to_store, false, parallel, grain_size, "", seed_in,
                  data_created_in_get_measurement_covariance_estimators);
              // (dev/plan.md 3.21 (c2''')) SMC^2 with an ensemble Kalman filter
              if (factor_plan.data_tempered)
              {
                inner_enkf->set_follows_caller_target(plan.data_variable, data_tempering.anneal);
                plan.followed_observations.push_back(inner_enkf->observations());
              }
              // (C267) the exchange step: SMC^2 only, as for a particle filter
              double exchange_below = enkf_params[5];
              if (exchange_below >= 0.0)
              {
                if (!factor_plan.data_tempered)
                  stop("enkf(..., increase = exchange(...)) is for SMC^2: an outer SMC whose sequence = data(...) takes in the "
                       "filter's observations (in MCMC it would change the target of a running chain).");
                double exchange_multiply = enkf_params[6];
                size_t exchange_max = enkf_params[7];
                if (exchange_multiply <= 1.0)
                  stop("exchange(..., multiply = m): m must be more than 1.");
                inner_enkf->set_exchange(exchange_below, exchange_multiply, exchange_max);
              }
              new_likelihood_estimator = inner_enkf;
            } else if (iequals(llhd_type, "pf")) {
              RecipeValue pf_params =
                  get_pf_info(model_parameters, llhd_details, llhd_type);

              size_t number_of_particles = pf_params[0];
              size_t iterations_to_store = pf_params[1];
              bool parallel = pf_params[2];
              size_t grain_size = pf_params[3];
              // an inner model that calls R cannot run in parallel threads
              if (parallel && llhd_model.contains("uses_r_functions") &&
                  (llhd_model["uses_r_functions"]).as_bool())
                Rcpp::stop("A particle filter with parallel = TRUE needs a model written in C++ (or with ilike:: built-in functions): its model has functions written in R, Python or Julia, or R expressions, which cannot be called from parallel threads.");

              ParticleFilter *inner_pf = get_particle_filter(
                  rng_in, data_in, llhd_model, model_parameters,
                  algorithm_parameters, number_of_particles,
                  iterations_to_store, false, parallel, grain_size, "", seed_in,
                  data_created_in_get_likelihood_estimators,
                  data_created_in_get_measurement_covariance_estimators);
              // nested: keep its stored generations compactly and write them
              // concatenated over the outer points
              inner_pf->set_compact_stored_output(true);
              // (dev/plan.md 3.21 (c2')) SMC^2: the caller's data tempering
              // takes in the filter's observations one at a time
              if (factor_plan.data_tempered) {
                // (a filter takes its observations in time order, whatever
                // data(..., order = ...) says: that is for factors whose
                // terms may enter in any order)
                // (3.21 (c2''), F1) with anneal = true the caller may take in
                // a share of an observation; the SMC builder makes each whole
                // observation a point of the schedule (followed_observations)
                inner_pf->set_follows_caller_target(plan.data_variable, data_tempering.anneal);
                plan.followed_observations.push_back(inner_pf->observations());
              }
              // (3.21 (c2''), F2) the exchange step: SMC^2 only (an outer SMC
              // whose data tempering takes in the filter's observations)
              double exchange_below = pf_params[4];
              if (exchange_below >= 0.0)
              {
                if (!factor_plan.data_tempered)
                  stop("pf(..., increase = exchange(...)) is for SMC^2: an outer SMC whose sequence = data(...) takes in the "
                       "filter's observations (in MCMC it would change the target of a running chain).");
                double exchange_multiply = pf_params[5];
                size_t exchange_max = pf_params[6];
                if (exchange_multiply <= 1.0)
                  stop("exchange(..., multiply = m): m must be more than 1.");
                inner_pf->set_exchange(exchange_below, exchange_multiply, exchange_max);
              }
              // threads for an inner particle filter with parallel = TRUE
              if (model.contains("number_of_threads"))
                inner_pf->set_number_of_threads(size_t((model["number_of_threads"]).as_double()));
              new_likelihood_estimator = inner_pf;
            } else if (iequals(llhd_type, "is")) {
              // integrate u with is(points = N) (dev/plan.md 3.21, step
              // (b)): an importance sampler over u, given the outer
              // parameters, on the inner model the lowering made (u's prior,
              // the proposal, and the factors of the data that depend on u)
              RecipeValue is_params = get_is_info(model_parameters, llhd_details);
              // (dev/plan.md 3.21, step (b)) when the outer sequence changes
              // this factor (an ABC tolerance read by the data's kernel), the
              // inner points are kept and weighted again under each new outer
              // target, as in the ABC-SMC of arXiv:2211.02172 (alg. 3)
              NestedSequence nested;
              if (!factor_is_fixed_in_smc) {
                nested.types = sequencer_types;
                nested.variables = sequencer_variables;
              }
              size_t number_of_points = is_params[0];
              bool parallel = is_params[1];
              size_t grain_size = is_params[2];
              new_likelihood_estimator = get_importance_sampler(
                  rng_in, data_in, llhd_model, model_parameters, algorithm_parameters,
                  number_of_points, parallel, grain_size, "", seed_in,
                  data_created_in_get_likelihood_estimators,
                  data_created_in_get_measurement_covariance_estimators, &nested);
            } else if (iequals(llhd_type, "smc")) {
              // integrate u with smc(particles = N, move = …, sequence = …)
              // (dev/plan.md 3.21, (c1)): an SMC sampler with MCMC moves over
              // u, given the outer parameters, through its own sequence (e.g.
              // decreasing ABC tolerances: Prangle et al.'s rare-event ABC);
              // its evidence estimate is the factor's likelihood. The methods
              // were prepared by the R binding (resolve_integrations).
              if (!llhd_details.contains("methods"))
                stop("integrate ... with smc(...): the inner algorithm was not prepared (compile the file with ilike::compile).");
              RecipeValue smc_params = get_integrate_smc_info(model_parameters, llhd_details);
              size_t number_of_particles = smc_params[0];
              bool parallel = smc_params[1];
              size_t grain_size = smc_params[2];
              RecipeValue methods = llhd_details["methods"];
              // (dev/plan.md 3.21 (c2)) an outer sequence that changes this
              // factor (a decreasing tolerance): the inner SMC follows it,
              // one step per outer step (rare-event ABC-SMC^2,
              // arXiv:2211.02172); otherwise it runs its own sequence
              // (rare-event ABC-MCMC, Prangle et al.)
              std::string u_name = llhd_details["model_file"].as_string();
              bool own_sequence = methods["smc_sequence"].contains("sequences");
              NestedSequence nested;
              if (!factor_is_fixed_in_smc) {
                if (own_sequence)
                  stop("integrate " + u_name + " with smc(...): the outer SMC's sequence changes this factor, so the "
                       "inner SMC follows the outer targets: remove sequence (and target) from smc(...).");
                nested.types = sequencer_types;
                nested.variables = sequencer_variables;
              } else if (!own_sequence) {
                stop("integrate " + u_name + " with smc(...): the SMC's sequence is needed (e.g. sequence = "
                     "parameter(epsilon, [Inf, 2, 1, 0.5])), unless an outer SMC sequence changes this factor.");
              }
              new_likelihood_estimator = get_smc_mcmc_move(
                  rng_in, data_in, llhd_model, model_parameters, algorithm_parameters,
                  number_of_particles, methods["mcmc_termination"], methods["mcmc_weights"],
                  methods["smc_sequence"], methods["adaptive_target"], methods["smc_termination"],
                  2, false, parallel, grain_size, "", seed_in,
                  data_created_in_get_likelihood_estimators,
                  data_created_in_get_measurement_covariance_estimators, &nested);
              // (dev/plan.md 3.21 (c4)) early stopping in an outer MCMC
              // move: ABC kernels (the lowering's flag), tolerances not increasing
              if (nested.types.empty())
              {
                static_cast<SMCMCMCMove*>(new_likelihood_estimator)->set_stops_early(
                    llhd_details.contains("kernels_only") && llhd_details["kernels_only"].as_bool());
                // (dev/plan.md 3.21 (c6), step 1) its own sequence (rare-event
                // ABC-MCMC): the stored generations kept compactly, as a nested
                // filter's (C100), so an outer chain's stored iterations do not
                // each hold a rich copy of the inner run. (An inner SMC that
                // follows the outer sequence keeps rich states between steps.)
                static_cast<SMCMCMCMove*>(new_likelihood_estimator)->set_compact_stored_output(true);
              }
            } else {
              stop("Invalid type for 'likelihood'.");
            }

            full_index_vector.push_back(likelihood_estimators.size());
            if (factor_plan.in_weight)
              without_cancelled_index_vector.push_back(
                  likelihood_estimators.size());
            if (factor_plan.in_without_priors)
              without_priors_index_vector.push_back(likelihood_estimators.size());

            if (factor_is_annealed) {
              likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                  rng_in, seed_in, data_in, new_likelihood_estimator, power,
                  annealing_variable, false));
            } else {
              likelihood_estimators.push_back(new_likelihood_estimator);
            }
          } else {
            stop("likelihood factor not well specified.");
          }
        } else if (current_factor.contains("unnormalised_normaliser")) {
          // (dev/plan.md 3.20, S4) an estimate of 1/Z(theta) of an
          // unnormalised distribution (normaliser = sav(...))
          RecipeValue entry = current_factor["unnormalised_normaliser"];
          std::vector<std::string> normaliser = entry["normaliser"];
          std::vector<double> at = entry["at"];
          UnnormalisedExactFunctionPtr exact = NULL;
          if (entry.contains("exact"))
            exact = entry["exact"].function_pointer<UnnormalisedExactFunctionPtr>();
          LikelihoodEstimator *new_likelihood_estimator = new UnnormalisedNormaliserEstimator(
              rng_in, seed_in, data_in, entry["variable"].as_string(),
              entry["potential"].function_pointer<UnnormalisedPotentialPtr>(), exact,
              entry["sweep"].function_pointer<UnnormalisedSweepFunctionPtr>(),
              size_t(entry["sweeps"].as_double()), size_t(entry["points"].as_double()),
              normaliser, at, factor_is_fixed_in_smc);
          full_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_weight)
            without_cancelled_index_vector.push_back(likelihood_estimators.size());
          if (factor_plan.in_without_priors)
            without_priors_index_vector.push_back(likelihood_estimators.size());
          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.contains("simulate_data_model")) {
          // sbi(simulate = …, data = …) without a method: only MFDS uses a
          // simulator directly
          stop("Factor " + std::to_string(i + 1) + " is a simulator without a method (sbi(simulate = …, data = …)): "
               "only MFDS can use it. For other algorithms, give the simulation-based method, e.g. "
               "sbi(sl(100), simulate = …, data = …).");
        } else {
          stop("Invalid factor.");
        }
      } else {
        stop("Error in factors section of model file.");
      }
    }
  }

  // Default annealing (all factors) with an importance proposal q: the target
  // path is (all factors)^p q^(1-p), starting at q. When only chosen factors
  // are annealed, the path is (other factors) (chosen factors)^p, with no q
  // term; the first weights are then (target / q) as usual.
  if (plan.add_proposal_path_factor) {
    full_index_vector.push_back(likelihood_estimators.size());
    without_cancelled_index_vector.push_back(likelihood_estimators.size());

    ExactLikelihoodEstimator *proposal_for_evaluation =
        new ExactLikelihoodEstimator(rng_in, seed_in, data_in, proposal_in,
                                     true);
    likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
        rng_in, seed_in, data_in, proposal_for_evaluation, second_power,
        annealing_variable, false));
  }

  full_index = new VectorIndex(full_index_vector);
  without_cancelled_index = new VectorIndex(without_cancelled_index_vector);
  without_priors_index = new VectorIndex(without_priors_index_vector);

  return likelihood_estimators;

}

std::vector<KalmanUpdater *> get_kalman_updaters(const RecipeValue &model,
                                                 const RecipeValue &model_parameters) {
  std::vector<KalmanUpdater *> updaters;

  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];

        if (is_linear_gaussian_data_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_data_model(current_factor, model_parameters);
          size_t type = lg.is_function ? 2 : 1;

          if (type == 1) {
            updaters.push_back(new ExactKalmanUpdater(
                lg.state_variable_string,
                lg.variable,
                lg.matrix,
                lg.covariance));
          } else if (type == 2) {
            updaters.push_back(new ExactKalmanUpdater(
                lg.state_variable_string,
                lg.variable,
                lg.matrix_function,
                lg.covariance_function));
          } else {
            stop("linear_gaussian specificiation is invalid.");
          }
        } else if (current_factor.contains(
                       "nonlinear_gaussian_data_covariance")) {
          stop("Unscented Kalman filter not yet supported.");

        }
      } else {
        stop("Error in factors section of model file.");
      }
    }
  } else {
    stop("No factors found in model file.");
  }

  if (updaters.size() == 0)
    stop("Valid Kalman updater not found.");
  else
    return updaters;
}

KalmanPredictor *get_kalman_predictor(const RecipeValue &model,
                                      const RecipeValue &model_parameters) {

  if (model.contains("transition_model")) {
    RecipeValue factors = model["transition_model"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];

        if (is_linear_gaussian_transition_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_transition_model(
              current_factor, model_parameters, false);
          if (!lg.is_function) {
            return new ExactKalmanPredictor(lg.matrix, lg.covariance);
          } else {
            return new ExactKalmanPredictor(lg.matrix_function,
                                            lg.covariance_function);
          }
        } else if (current_factor.contains(
                       "nonlinear_gaussian_transition_covariance")) {
          stop("Unscented Kalman filter not yet supported.");

        }
      } else {
        stop("Error in transition_models section of model file.");
      }
    }
  } else {
    stop("No transition_models found in model file.");
  }

  stop("Valid Kalman predictor not found.");
}

// (C269, PMMH P1) a transition's reads = [...] (what its matrix and covariance
// functions read): never the state variable nor the time index's variables
static std::vector<std::string> transition_declared_reads(const RecipeValue &model, const RecipeValue &model_parameters,
                                                          const RecipeValue &factor, const std::string &key,
                                                          const std::string &state_variable)
{
  std::vector<std::string> reads;
  if (!factor.contains(key))
    return reads;
  reads = factor[key].as_strings();
  std::set<std::string> forbidden;
  forbidden.insert(state_variable);
  if (model.contains("method"))
  {
    RecipeValue methods = model["method"];
    for (size_t i = 0; i < methods.size(); ++i)
      if (methods[i].contains("ssm"))
      {
        RecipeValue ssm = get_ssm_info(model, model_parameters);
        forbidden.insert(ssm[0].as_string()); // (the index)
        forbidden.insert(ssm[3].as_string()); // (the time)
        forbidden.insert(ssm[5].as_string()); // (the time difference)
        break;
      }
  }
  for (const std::string &name : reads)
    if (forbidden.count(name))
      stop("a transition's reads = [...] may name only parameters, not the state (" + state_variable +
           ") or the time index's variables: " + name + " is one of them (its value changes within a filter run).");
  return reads;
}

ProposalKernel *get_transition_proposal(const RecipeValue &model,
                                        const RecipeValue &model_parameters) {

  if (model.contains("transition_proposal")) {
    RecipeValue factors = model["transition_proposal"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];

        if (current_factor.contains(
                "linear_gaussian_transition_covariance")) {
          if (!current_factor.contains(
                  "linear_gaussian_transition_matrix")) {
            stop("linear_gaussian must contain both matrix and covariance.");
          }

          if (!current_factor.contains("type")) {
            stop("linear_gaussian must contain a type.");
          }

          if (!current_factor.contains(
                  "linear_gaussian_transition_variable")) {
            stop("linear_gaussian must contain "
                 "linear_gaussian_transition_variable.");
          }

          RecipeValue linear_gaussian_transition_matrix_value =
              current_factor["linear_gaussian_transition_matrix"];
          RecipeValue linear_gaussian_transition_covariance_value =
              current_factor["linear_gaussian_transition_covariance"];
          RecipeValue linear_gaussian_transition_variable_value =
              current_factor["linear_gaussian_transition_variable"];

          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new LinearGaussianNoiseProposalKernel(
                (linear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                (linear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                (linear_gaussian_transition_matrix_value).function_pointer<GetNoParamsMatrixPtr>()(),
                (linear_gaussian_transition_covariance_value).function_pointer<GetNoParamsMatrixPtr>()());
          } else if (type == 2) {
            LinearGaussianNoiseFunctionProposalKernel *kernel = new LinearGaussianNoiseFunctionProposalKernel(
                (linear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                (linear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                (linear_gaussian_transition_matrix_value).function_pointer<GetMatrixPtr>(),
                (linear_gaussian_transition_covariance_value).function_pointer<GetMatrixPtr>());
            // (C269, PMMH P1) a constant covariance beside a matrix function
            if (current_factor.contains("linear_gaussian_transition_covariance_constant") &&
                current_factor["linear_gaussian_transition_covariance_constant"].as_bool())
              kernel->set_covariance_constant();
            kernel->set_declared_reads(transition_declared_reads(model, model_parameters, current_factor,
                                                                 "linear_gaussian_transition_reads", (linear_gaussian_transition_variable_value).function_pointer<StringPtr>()()));
            proposal = kernel;
          } else {
            stop("linear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.contains(
                       "nonlinear_gaussian_transition_covariance")) {

          if (!current_factor.contains(
                  "nonlinear_gaussian_transition_function")) {
            stop("nonlinear_gaussian_transition must contain both function and "
                 "covariance.");
          }

          if (!current_factor.contains("type")) {
            stop("nonlinear_gaussian_transition must contain a type.");
          }

          if (!current_factor.contains(
                  "nonlinear_gaussian_transition_variable")) {
            stop("nonlinear_gaussian_transition must contain "
                 "nonlinear_gaussian_transition_variable.");
          }

          RecipeValue nonlinear_gaussian_transition_function_value =
              current_factor["nonlinear_gaussian_transition_function"];
          RecipeValue nonlinear_gaussian_transition_covariance_value =
              current_factor["nonlinear_gaussian_transition_covariance"];
          RecipeValue nonlinear_gaussian_transition_variable_value =
              current_factor["nonlinear_gaussian_transition_variable"];
          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>((nonlinear_gaussian_transition_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_transition_covariance_value).function_pointer<GetNoParamsMatrixPtr>()());
          } else if (type == 2) {
            NonLinearGaussianNoiseFunctionProposalKernel *kernel = new NonLinearGaussianNoiseFunctionProposalKernel(
                (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>((nonlinear_gaussian_transition_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_transition_covariance_value).function_pointer<GetMatrixPtr>());
            // (C269, PMMH P1) a constant covariance beside a function
            if (current_factor.contains("nonlinear_gaussian_transition_covariance_constant") &&
                current_factor["nonlinear_gaussian_transition_covariance_constant"].as_bool())
              kernel->set_covariance_constant();
            kernel->set_declared_reads(transition_declared_reads(model, model_parameters, current_factor,
                                                                 "nonlinear_gaussian_transition_reads", (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()()));
            proposal = kernel;
          } else {
            stop("nonlinear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.contains("transition_proposal")) {

          if ((current_factor["transition_proposal"]).is_list()) {

            RecipeValue current_transition_proposal =
                current_factor["transition_proposal"];
            if (current_transition_proposal.contains("type") &&
                current_transition_proposal.contains("variables")) {
              return make_builtin_symmetric_proposal(current_transition_proposal, BuiltinRole::TransitionProposal, model_parameters);
            } else {
              stop("Missing information in transition model section of model "
                   "file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.contains(
                       "simulate_transition_proposal")) {
          RecipeValue simulate_transition_proposal_value =
              current_factor["simulate_transition_proposal"];

          if (current_factor.contains(
                  "evaluate_log_transition_proposal")) {
            RecipeValue evaluate_log_transition_proposal_value =
                current_factor["evaluate_log_transition_proposal"];
            return new CustomGuidedNoParamsProposalKernel(
                (simulate_transition_proposal_value).function_pointer<SimulateGuidedNoParamsMCMCProposalPtr>(),
                (evaluate_log_transition_proposal_value).function_pointer<EvaluateLogGuidedNoParamsMCMCProposalPtr>());
          } else {
            return new CustomGuidedNoParamsProposalKernel(
                (simulate_transition_proposal_value).function_pointer<SimulateGuidedNoParamsMCMCProposalPtr>());
          }
        }
      } else {
        stop("Error in transition_proposals section of model file.");
      }
    }
  } else {
    return NULL;
  }

  stop("No valid transition proposal found.");
}

ProposalKernel *get_transition_model(const RecipeValue &model,
                                     const RecipeValue &model_parameters) {

  if (model.contains("transition_model")) {
    RecipeValue factors = model["transition_model"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];

        if (is_linear_gaussian_transition_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_transition_model(
              current_factor, model_parameters, true);
          ProposalKernel *proposal;
          if (!lg.is_function) {
            proposal = new LinearGaussianNoiseProposalKernel(
                lg.variable, lg.variable, lg.matrix, lg.covariance);
          } else {
            LinearGaussianNoiseFunctionProposalKernel *kernel = new LinearGaussianNoiseFunctionProposalKernel(
                lg.variable, lg.variable, lg.matrix_function,
                lg.covariance_function);
            // (C269, PMMH P1) a constant covariance beside a matrix function
            if (current_factor.contains("linear_gaussian_transition_covariance_constant") &&
                current_factor["linear_gaussian_transition_covariance_constant"].as_bool())
              kernel->set_covariance_constant();
            kernel->set_declared_reads(transition_declared_reads(model, model_parameters, current_factor,
                                                                 "linear_gaussian_transition_reads", lg.variable));
            proposal = kernel;
          }

          return proposal;
        } else if (current_factor.contains(
                       "nonlinear_gaussian_transition_covariance")) {

          if (!current_factor.contains(
                  "nonlinear_gaussian_transition_function")) {
            stop("nonlinear_gaussian_transition must contain both function and "
                 "covariance.");
          }

          if (!current_factor.contains("type")) {
            stop("nonlinear_gaussian_transition must contain a type.");
          }

          if (!current_factor.contains(
                  "nonlinear_gaussian_transition_variable")) {
            stop("nonlinear_gaussian_transition must contain "
                 "nonlinear_gaussian_transition_variable.");
          }

          RecipeValue nonlinear_gaussian_transition_function_value =
              current_factor["nonlinear_gaussian_transition_function"];
          RecipeValue nonlinear_gaussian_transition_covariance_value =
              current_factor["nonlinear_gaussian_transition_covariance"];
          RecipeValue nonlinear_gaussian_transition_variable_value =
              current_factor["nonlinear_gaussian_transition_variable"];
          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>((nonlinear_gaussian_transition_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_transition_covariance_value).function_pointer<GetNoParamsMatrixPtr>()());
          } else if (type == 2) {
            NonLinearGaussianNoiseFunctionProposalKernel *kernel = new NonLinearGaussianNoiseFunctionProposalKernel(
                (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()(),
                std::make_shared<Transform>((nonlinear_gaussian_transition_function_value).function_pointer<TransformPtr>()),
                (nonlinear_gaussian_transition_covariance_value).function_pointer<GetMatrixPtr>());
            // (C269, PMMH P1) a constant covariance beside a function
            if (current_factor.contains("nonlinear_gaussian_transition_covariance_constant") &&
                current_factor["nonlinear_gaussian_transition_covariance_constant"].as_bool())
              kernel->set_covariance_constant();
            kernel->set_declared_reads(transition_declared_reads(model, model_parameters, current_factor,
                                                                 "nonlinear_gaussian_transition_reads", (nonlinear_gaussian_transition_variable_value).function_pointer<StringPtr>()()));
            proposal = kernel;
          } else {
            stop("nonlinear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.contains("transition_model")) {
          if ((current_factor["transition_model"]).is_list()) {

            RecipeValue current_transition_model = current_factor["transition_model"];
            if (current_transition_model.contains("type") &&
                current_transition_model.contains("variables")) {
              return make_builtin_symmetric_proposal(current_transition_model, BuiltinRole::TransitionModel, model_parameters);
            } else {
              stop("Missing information in transition model section of model "
                   "file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.contains(
                       "simulate_transition_model")) {
          RecipeValue simulate_transition_model_value =
              current_factor["simulate_transition_model"];

          if (current_factor.contains(
                  "evaluate_log_transition_model")) {
            RecipeValue evaluate_log_transition_model_value =
                current_factor["evaluate_log_transition_model"];
            return new CustomNoParamsProposalKernel(
                (simulate_transition_model_value).function_pointer<SimulateNoParamsMCMCProposalPtr>(),
                (evaluate_log_transition_model_value).function_pointer<EvaluateLogNoParamsMCMCProposalPtr>());
          } else {
            return new CustomNoParamsProposalKernel(
                (simulate_transition_model_value).function_pointer<SimulateNoParamsMCMCProposalPtr>());
          }
        }
      } else {
        stop("Error in transition_models section of model file.");
      }
    }
  } else {
    stop("No transition_models found in model file.");
  }

  stop("No valid transition model found.");
}

void get_gaussian_prior_moments(const RecipeValue &model, const RecipeValue &model_parameters,
                                Parameters &means, Parameters &covariances,
                                std::vector<std::string> &log_variables) {
  const std::string requirement =
      "MFDS requires a Gaussian prior on every parameter, or a log-normal one (the exp of a Gaussian): "
      "theta ~ norm(mean, sd), mvnorm(mean, covariance), lnorm(meanlog, sdlog) or "
      "mvlnorm(meanlog, covlog), with arguments that do not depend on parameters.";
  if (!model.contains("factor"))
    stop(requirement + " No prior was found.");
  RecipeValue factors = model["factor"];
  size_t number_found = 0;
  for (size_t i = 0; i < factors.size(); ++i) {
    RecipeValue current_factor = factors[i];
    if (current_factor.contains("evaluate_log_prior") ||
        current_factor.contains("simulate_prior"))
      stop(requirement + " Factor " + std::to_string(i + 1) + " is a user-defined prior.");
    if (!current_factor.contains("prior"))
      continue;
    RecipeValue current_prior = current_factor["prior"];
    std::string type = current_prior["type"];
    BuiltinArguments arguments = bind_builtin(current_prior, BuiltinRole::Prior, model_parameters);
    std::string variable = arguments.variable("variable");
    if (iequals(type, "norm")) {
      double sd = arguments.real("sd");
      means[variable] = arma::mat(1, 1, arma::fill::value(arguments.real("mean")));
      covariances[variable] = arma::mat(1, 1, arma::fill::value(sd * sd));
    } else if (iequals(type, "mvnorm")) {
      means[variable] = arma::mat(arguments.vector("mean"));
      covariances[variable] = arguments.matrix("cov");
    } else if (iequals(type, "lnorm")) {
      // (the Gaussian is on log(variable))
      double sdlog = arguments.real("sdlog");
      means[variable] = arma::mat(1, 1, arma::fill::value(arguments.real("meanlog")));
      covariances[variable] = arma::mat(1, 1, arma::fill::value(sdlog * sdlog));
      log_variables.push_back(variable);
    } else if (iequals(type, "mvlnorm")) {
      means[variable] = arma::mat(arguments.vector("meanlog"));
      covariances[variable] = arguments.matrix("covlog");
      log_variables.push_back(variable);
    } else {
      stop(requirement + " Factor " + std::to_string(i + 1) +
           " is a built-in prior of type " + type + ".");
    }
    number_found = number_found + 1;
  }
  if (number_found == 0)
    stop(requirement + " No prior was found.");
}

RecipeValue get_prior_mean_and_covariance(const RecipeValue &model,
                                   const RecipeValue &model_parameters) {
  if (!model.contains("factor"))
    stop(kf_prior_requirement + " No prior was found.");

  RecipeValue factors = model["factor"];

  for (size_t i = 0; i < factors.size(); ++i) {

    if (!(factors[i]).is_list())
      stop("Error in factors section of model file.");

    RecipeValue current_factor = factors[i];

    if (current_factor.contains("evaluate_log_prior") ||
        current_factor.contains("simulate_prior"))
      stop(kf_prior_requirement + " Factor " + std::to_string(i + 1) +
           " is a user-defined prior.");

    if (!current_factor.contains("prior"))
      continue;

    if (!(current_factor["prior"]).is_list())
      stop("Error in factors section of model file.");

    RecipeValue current_prior = current_factor["prior"];
    if (!(current_prior.contains("type") &&
          current_prior.contains("variables")))
      stop("Missing information for prior in model file.");

    std::string type = current_prior["type"];

    if (iequals(type, "norm")) {
      BuiltinArguments arguments =
          bind_builtin(current_prior, BuiltinRole::Prior, model_parameters);
      double sd = arguments.real("sd");

      arma::colvec mean_vec(1);
      mean_vec[0] = arguments.real("mean");

      arma::mat cov(1, 1);
      cov(0, 0) = sd * sd;

      return RecipeValue::list({mean_vec, cov});
    } else if (iequals(type, "mvnorm")) {
      BuiltinArguments arguments =
          bind_builtin(current_prior, BuiltinRole::Prior, model_parameters);
      return RecipeValue::list({arguments.vector("mean"), arguments.matrix("cov")});
    } else {
      stop(kf_prior_requirement + " Factor " + std::to_string(i + 1) +
           " is a built-in prior of type " + type + ".");
    }
  }

  stop(kf_prior_requirement + " No prior was found.");
}

std::vector<MeasurementCovarianceEstimator *>
get_measurement_covariance_estimators(
    RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in,
    const RecipeValue &model, const RecipeValue &model_parameters,
    const std::vector<size_t> &factor_indices,
    std::shared_ptr<Transform> transform,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {

  std::vector<MeasurementCovarianceEstimator *>
      measurement_covariance_estimators;

  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];

    std::vector<size_t> indices_to_use;
    if (factor_indices.size() == 0) {
      indices_to_use.reserve(factors.length());
      for (size_t i = 0; i < factors.length(); ++i) {
        indices_to_use.push_back(i);
      }
    } else {
      indices_to_use = factor_indices;
    }

    for (size_t i = 0; i < indices_to_use.size(); ++i) {
      if (indices_to_use[i] >= factors.length()) {
        stop("Error in enk_likelihood_index: index does not index a specified "
             "factor.");
      }

      if ((factors[indices_to_use[i]]).is_list()) {
        RecipeValue current_factor = factors[indices_to_use[i]];
        SimulateModelPtr simulate_data_model_in;
        if (current_factor.contains("simulate_data_model")) {
          RecipeValue simulate_data_model_value = current_factor["simulate_data_model"];
          simulate_data_model_in =
              (simulate_data_model_value).function_pointer<SimulateModelPtr>();

          if (!current_factor.contains("data_variable")) {
            stop("simulator-based model must be associated with a "
                 "data_variable.");
          }

          RecipeValue data_variable_value = current_factor["data_variable"];

          std::string augmented_data_variable_names =
              (data_variable_value).function_pointer<StringPtr>()();
          std::vector<std::string> data_variable_names =
              split(augmented_data_variable_names, ';');

          MeasurementCovarianceEstimator *new_measurement_covariance_estimator;

          if (current_factor.contains("summary_statistics")) {
            RecipeValue summary_statistics_value = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

            data_created_in_get_measurement_covariance_estimators.push_back(
                summary_stats(*data_in));

            new_measurement_covariance_estimator =
                new GenericMeasurementCovarianceEstimator(
                    rng_in, seed_in,
                    &data_created_in_get_measurement_covariance_estimators
                         .back(),
                    transform, std::make_shared<Transform>(summary_stats),
                    simulate_data_model_in, data_variable_names);
          } else {

            new_measurement_covariance_estimator =
                new GenericMeasurementCovarianceEstimator(
                    rng_in, seed_in, data_in, transform, NULL,
                    simulate_data_model_in, data_variable_names);
          }

          measurement_covariance_estimators.push_back(
              new_measurement_covariance_estimator);
        } else if (is_linear_gaussian_data_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_data_model(current_factor, model_parameters);
          size_t type = lg.is_function ? 2 : 1;

          std::string augmented_state_variable_names =
              lg.state_variable_string;
          std::vector<std::string> state_variable_names =
              split(augmented_state_variable_names, ';');

          MeasurementCovarianceEstimator *new_measurement_covariance_estimator;

          if (current_factor.contains("summary_statistics")) {
            RecipeValue summary_statistics_value = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

            data_created_in_get_measurement_covariance_estimators.push_back(
                summary_stats(*data_in));

            if (type == 1) {
              new_measurement_covariance_estimator =
                  new DirectLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in,
                      &data_created_in_get_measurement_covariance_estimators
                           .back(),
                      transform, std::make_shared<Transform>(summary_stats),
                      lg.matrix,
                      lg.covariance,
                      lg.variable,
                      state_variable_names);
            } else if (type == 2) {
              new_measurement_covariance_estimator =
                  new DirectLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in,
                      //&data_created_in_get_measurement_covariance_estimators.back(),
                      transform,
                      // std::make_shared<Transform>(summary_stats),
                      NULL,
                      lg.matrix_function,
                      lg.covariance_function,
                      lg.variable,
                      state_variable_names);
            } else {
              stop("linear_gaussian_data specificiation is invalid.");
            }
          } else {

            if (type == 1) {
              new_measurement_covariance_estimator =
                  new DirectLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in, transform, NULL,
                      lg.matrix,
                      lg.covariance,
                      lg.variable,
                      state_variable_names);
            } else if (type == 2) {
              new_measurement_covariance_estimator =
                  new DirectLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in, transform, NULL,
                      lg.matrix_function,
                      lg.covariance_function,
                      lg.variable,
                      state_variable_names);
            } else {
              stop("linear_gaussian_data specificiation is invalid.");
            }
          }

          measurement_covariance_estimators.push_back(
              new_measurement_covariance_estimator);
        } else if (current_factor.contains(
                       "nonlinear_gaussian_data_covariance")) {
          if (!current_factor.contains(
                  "nonlinear_gaussian_data_function")) {
            stop("nonlinear_gaussian_data must contain both function and "
                 "covariance.");
          }

          if (!current_factor.contains("type")) {
            stop("nonlinear_gaussian_data must contain a type.");
          }

          RecipeValue nonlinear_gaussian_data_function_value =
              current_factor["nonlinear_gaussian_data_function"];
          RecipeValue nonlinear_gaussian_data_covariance_value =
              current_factor["nonlinear_gaussian_data_covariance"];
          RecipeValue nonlinear_gaussian_data_variable_value =
              current_factor["nonlinear_gaussian_data_variable"];
          size_t type = current_factor["type"];

          TransformPtr nonlinear_gaussian_data_function =
              (nonlinear_gaussian_data_function_value).function_pointer<TransformPtr>();

          MeasurementCovarianceEstimator *new_measurement_covariance_estimator;

          if (current_factor.contains("summary_statistics")) {
            RecipeValue summary_statistics_value = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                (summary_statistics_value).function_pointer<SummaryStatisticsPtr>();

            data_created_in_get_measurement_covariance_estimators.push_back(
                summary_stats(*data_in));

            if (type == 1) {
              new_measurement_covariance_estimator =
                  new DirectNonLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in,
                      &data_created_in_get_measurement_covariance_estimators
                           .back(),
                      transform, std::make_shared<Transform>(summary_stats),
                      std::make_shared<Transform>(
                          nonlinear_gaussian_data_function),
                      (nonlinear_gaussian_data_covariance_value).function_pointer<GetNoParamsMatrixPtr>()(),
                      (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()());
            } else if (type == 2) {
              new_measurement_covariance_estimator =
                  new DirectNonLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in,
                      //&data_created_in_get_measurement_covariance_estimators.back(),
                      transform,
                      // std::make_shared<Transform>(summary_stats),
                      NULL,
                      std::make_shared<Transform>(
                          nonlinear_gaussian_data_function),
                      (nonlinear_gaussian_data_covariance_value).function_pointer<GetMatrixPtr>(),
                      (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()());
            } else {
              stop("nonlinear_gaussian_data specificiation is invalid.");
            }
          } else {

            if (type == 1) {
              new_measurement_covariance_estimator =
                  new DirectNonLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in, transform, NULL,
                      std::make_shared<Transform>(
                          nonlinear_gaussian_data_function),
                      (nonlinear_gaussian_data_covariance_value).function_pointer<GetNoParamsMatrixPtr>()(),
                      (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()());
            } else if (type == 2) {
              new_measurement_covariance_estimator =
                  new DirectNonLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in, transform, NULL,
                      std::make_shared<Transform>(
                          nonlinear_gaussian_data_function),
                      (nonlinear_gaussian_data_covariance_value).function_pointer<GetMatrixPtr>(),
                      (nonlinear_gaussian_data_variable_value).function_pointer<StringPtr>()());
            } else {
              stop("nonlinear_gaussian_data specificiation is invalid.");
            }
          }

          measurement_covariance_estimators.push_back(
              new_measurement_covariance_estimator);
        }
      }
      // else
      //{
      //   stop("For use in an ensemble Kalman method, the factor must be either
      //   linear Gaussian, non-linear Gaussian, or simulate_model.");
      // }
    }
  }

  return measurement_covariance_estimators;

}
