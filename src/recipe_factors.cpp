#include "recipe_internal.h"

// Error for annealing a prior that is also the importance proposal: it
// cancels in the weights, so it is not part of the target and its power would
// be ignored.
std::string prior_cannot_be_annealed(size_t factor_index) {
  return "smc_sequence annealing: factor " + std::to_string(factor_index + 1) +
         " is a prior that is also used as the importance proposal, so it "
         "cancels and cannot be annealed. Remove it from the factors to "
         "anneal, or give an is_proposal.";
}

std::vector<LikelihoodEstimator *> get_likelihood_estimators(
    RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in,
    const List &model, const List &model_parameters,
    const List &algorithm_parameter_list, bool include_priors,
    const std::vector<std::string> &sequencer_types,
    const std::vector<std::string> &sequencer_variables,
    const std::vector<std::vector<double>> &sequencer_schedules,
    IndependentProposalKernel *proposal_in,
    VectorIndex *&without_cancelled_index, VectorIndex *&without_priors_index,
    VectorIndex *&full_index, bool &any_annealing,
    const std::vector<int> &factors_affected_by_smc_sequence,
    std::vector<Data> &data_created_in_get_likelihood_estimators,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {

  std::vector<size_t> without_cancelled_index_vector;
  std::vector<size_t> without_priors_index_vector;
  std::vector<size_t> full_index_vector;

  any_annealing = false;
  std::string annealing_variable;
  for (size_t i = 0; i < sequencer_types.size(); ++i) {

    // if ( ( (sequencer_types[i]=="annealing") ||
    // (sequencer_types[i]=="tempering") ) && (!(
    // (sequencer_schedules[i].size()==2) && (sequencer_schedules[i][0]==0.0) &&
    // (sequencer_schedules[i][1]==1.0)) ) )
    if ((iequals(sequencer_types[i], "annealing")) ||
        (iequals(sequencer_types[i], "tempering"))) {
      if (any_annealing == false) {
        any_annealing = true;
        annealing_variable = sequencer_variables[i];
        break;
      } else {
        stop("Two annealing variables specified: not currently supported.");
      }
    }
  }

  PowerFunctionPtr power = annealing_power;
  PowerFunctionPtr second_power = annealing_one_minus_power;

  std::vector<LikelihoodEstimator *> likelihood_estimators;

  if (model.containsElementNamed("factor")) {

    List factors = model["factor"];

    for (size_t j = 0; j < factors_affected_by_smc_sequence.size(); ++j) {
      if ((factors_affected_by_smc_sequence[j] < 0) ||
          (size_t(factors_affected_by_smc_sequence[j]) >= factors.size()))
        stop("smc_sequence: the factors affected must be between 1 and " +
             std::to_string(factors.size()) + ".");
    }

    for (size_t i = 0; i < factors.size(); ++i) {

      bool factor_is_listed = false;
      for (size_t j = 0; j < factors_affected_by_smc_sequence.size(); ++j) {
        if (int(i) == factors_affected_by_smc_sequence[j])
          factor_is_listed = true;
      }
      // Annealing: the listed factors are raised to the power (all factors if
      // none are listed); the others enter the target unannealed. The value of
      // an annealed factor does not depend on the power, so it is evaluated
      // once (smcfixed) and only the power changes.
      // Other sequences (ilike::parameter): the listed factors are the ones
      // that depend on the sequence variable, so they are re-evaluated when
      // the target changes; with no list, every factor is. Without an SMC
      // sequence, factors are fixed.
      bool factor_is_annealed =
          any_annealing &&
          (factors_affected_by_smc_sequence.empty() || factor_is_listed);
      bool factor_changes_with_sequence =
          (!sequencer_types.empty()) &&
          (factors_affected_by_smc_sequence.empty() || factor_is_listed);
      bool factor_is_fixed_in_smc =
          any_annealing || !factor_changes_with_sequence;

      if (Rf_isNewList(factors[i])) {

        List current_factor = factors[i];
        if (current_factor.containsElementNamed("prior")) {
          DistributionFactor *new_factor;

          if (Rf_isNewList(current_factor["prior"])) {

            List current_prior = current_factor["prior"];
            if (current_prior.containsElementNamed("type") &&
                current_prior.containsElementNamed("variables")) {
              std::string type = current_prior["type"];

              if (iequals(type, "norm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                new_factor = new GaussianDistributionFactor(variable, mean, sd);
              } else if (iequals(type, "mvnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                new_factor =
                    new GaussianDistributionFactor(variable, mean, cov);
              } else if (iequals(type, "lnorm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                new_factor =
                    new LogGaussianDistributionFactor(variable, mean, sd);
              } else if (iequals(type, "mvlnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                new_factor =
                    new LogGaussianDistributionFactor(variable, mean, cov);
              } else if (iequals(type, "gamma")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];

                new_factor = new GammaDistributionFactor(variable, shape, rate);
              } else {
                Rcout << "Prior type " << type;
                stop("Prior type unknown.");
              }

              LikelihoodEstimator *new_likelihood_estimator =
                  new ExactLikelihoodEstimator(rng_in, seed_in, data_in,
                                               new_factor, true);

              full_index_vector.push_back(likelihood_estimators.size());
              if (include_priors) {
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
                if (any_annealing && factor_is_listed)
                  stop(prior_cannot_be_annealed(i));
                likelihood_estimators.push_back(new_likelihood_estimator);
              }
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.containsElementNamed("evaluate_log_prior")) {
          SEXP evaluate_log_prior_SEXP = current_factor["evaluate_log_prior"];
          // The gradient is optional: it is only needed by gradient-based
          // proposals (e.g. Langevin, Barker, HMC).
          DistributionFactor *new_factor;
          if (current_factor.containsElementNamed(
                  "evaluate_gradient_log_prior")) {
            SEXP evaluate_gradient_log_prior_SEXP =
                current_factor["evaluate_gradient_log_prior"];
            new_factor = new CustomDistributionFactor(
                load_evaluate_log_distribution(evaluate_log_prior_SEXP),
                load_evaluate_gradient_log_distribution(
                    evaluate_gradient_log_prior_SEXP));
          } else {
            new_factor = new CustomDistributionFactor(
                load_evaluate_log_distribution(evaluate_log_prior_SEXP));
          }
          LikelihoodEstimator *new_likelihood_estimator =
              new ExactLikelihoodEstimator(rng_in, seed_in, data_in, new_factor,
                                           factor_is_fixed_in_smc);
          full_index_vector.push_back(likelihood_estimators.size());

          if (include_priors) {

            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());

            if (factor_is_annealed) {

              likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                  rng_in, seed_in, data_in, new_likelihood_estimator, power,
                  annealing_variable, false));
            } else {

              likelihood_estimators.push_back(new_likelihood_estimator);
            }
          } else {
            if (any_annealing && factor_is_listed)
              stop(prior_cannot_be_annealed(i));
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.containsElementNamed("simulate_prior")) {
        } else if (current_factor.containsElementNamed(
                       "evaluate_log_likelihood")) {

          SEXP evaluate_log_likelihood_SEXP =
              current_factor["evaluate_log_likelihood"];
          // The gradient is optional: it is only needed by gradient-based
          // proposals (e.g. Langevin, Barker, HMC).
          LikelihoodFactor *new_factor;
          if (current_factor.containsElementNamed(
                  "evaluate_gradient_log_likelihood")) {
            SEXP evaluate_gradient_log_likelihood_SEXP =
                current_factor["evaluate_gradient_log_likelihood"];
            new_factor = new CustomLikelihoodFactor(
                load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP),
                load_evaluate_gradient_log_likelihood(
                    evaluate_gradient_log_likelihood_SEXP),
                data_in);
          } else {
            new_factor = new CustomLikelihoodFactor(
                load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP),
                data_in);
          }
          LikelihoodEstimator *new_likelihood_estimator =
              new ExactLikelihoodEstimator(rng_in, seed_in, data_in, new_factor,
                                           factor_is_fixed_in_smc);

          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(
              likelihood_estimators.size());
          without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.containsElementNamed("sbi_likelihood")) {
          SimulateModelPtr simulate_data_model_in;
          if (current_factor.containsElementNamed("simulate_data_model")) {
            SEXP simulate_data_model_SEXP =
                current_factor["simulate_data_model"];
            simulate_data_model_in =
                load_simulate_data_model(simulate_data_model_SEXP);
          } else {
            stop(
                "sbi_likelihood factor must also contain simulate_data_model.");
          }

          std::vector<std::string> data_variables_in;
          if (current_factor.containsElementNamed("data_variable")) {
            SEXP data_variable_SEXP = current_factor["data_variable"];
            std::string augmented_data_variable_names =
                load_string(data_variable_SEXP);
            data_variables_in = split(augmented_data_variable_names, ';');
          } else {
            stop("sbi_likelihood factor must also contain data_variable.");
          }

          LikelihoodEstimator *new_likelihood_estimator;

          List current_sbi = current_factor["sbi_likelihood"];
          if (current_sbi.containsElementNamed("type") &&
              current_sbi.containsElementNamed("parameters")) {
            std::string type = current_sbi["type"];

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
                List info = get_abc_euclidean_uniform_parameter_info(
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

                parallel = info[6];
                grain_size = info[7];
              } else {
                List info = get_abc_lp_uniform_parameter_info(
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

              if (current_factor.containsElementNamed("summary_statistics")) {
                SEXP summary_statistics_SEXP =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    load_summary_statistics(summary_statistics_SEXP);

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
              List info = get_abc_euclidean_uniform_parameter_info(
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

              parallel = info[6];
              grain_size = info[7];

              bool adaptive = false;
              for (auto k = sequencer_variables.begin();
                   k != sequencer_variables.end(); ++k) {
                if (*k == tolerance_variable) {
                  adaptive = true;
                }
              }

              if (current_factor.containsElementNamed("summary_statistics")) {
                SEXP summary_statistics_SEXP =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    load_summary_statistics(summary_statistics_SEXP);

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

              List info = get_abc_enki_parameter_info(model_parameters,
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

              EnsembleShifter *shifter;
              if (iequals(shifter_name, "stochastic")) {
                shifter = new StochasticEnsembleShifter();
              } else if (iequals(shifter_name, "sqrt")) {
                shifter = new SquareRootEnsembleShifter();
              } else if (iequals(shifter_name, "adjustment")) {
                shifter = new AdjustmentEnsembleShifter();
              } else {
                Rcout << "Shifter type: " << shifter_name << std::endl;
                Rcpp::stop("Invalid shifter type for EnKI");
              }

              bool adaptive = false;
              for (auto k = sequencer_variables.begin();
                   k != sequencer_variables.end(); ++k) {
                if (*k == tolerance_variable) {
                  adaptive = true;
                }
              }

              if (current_factor.containsElementNamed("summary_statistics")) {
                SEXP summary_statistics_SEXP =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    load_summary_statistics(summary_statistics_SEXP);

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
              List info =
                  get_sl_parameter_info(model_parameters, current_sbi, type);

              number_of_points = info[0];
              unbiased = info[1];
              store_output = info[2];
              parallel = info[3];
              grain_size = info[4];

              if (current_factor.containsElementNamed("summary_statistics")) {
                SEXP summary_statistics_SEXP =
                    current_factor["summary_statistics"];
                SummaryStatisticsPtr summary_stats =
                    load_summary_statistics(summary_statistics_SEXP);

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
          without_cancelled_index_vector.push_back(
              likelihood_estimators.size());
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
          without_cancelled_index_vector.push_back(
              likelihood_estimators.size());
          without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.containsElementNamed(
                       "nonlinear_gaussian_data_covariance")) {
          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_data_function")) {
            stop("nonlinear_gaussian_data must contain both function and "
                 "covariance.");
          }

          if (!current_factor.containsElementNamed("type")) {
            stop("nonlinear_gaussian_data must contain a type.");
          }

          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_data_variable")) {
            stop("nonlinear_gaussian_data must contain "
                 "nonlinear_gaussian_data_variable.");
          }

          // std::vector<std::string> variables =
          // data_in->get_vector_variables(); if (variables.size()!=1)
          //   stop("For this linear-Gaussian likelihood we can only have one
          //   variable present in the data.");

          SEXP nonlinear_gaussian_data_function_SEXP =
              current_factor["nonlinear_gaussian_data_function"];
          SEXP nonlinear_gaussian_data_covariance_SEXP =
              current_factor["nonlinear_gaussian_data_covariance"];

          SEXP nonlinear_gaussian_data_variable_SEXP =
              current_factor["nonlinear_gaussian_data_variable"];

          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                load_string(nonlinear_gaussian_data_variable_SEXP),
                std::make_shared<Transform>(
                    load_transform(nonlinear_gaussian_data_function_SEXP)),
                load_matrix(nonlinear_gaussian_data_covariance_SEXP));
          } else if (type == 2) {
            proposal = new NonLinearGaussianNoiseFunctionProposalKernel(
                load_string(nonlinear_gaussian_data_variable_SEXP),
                std::make_shared<Transform>(
                    load_transform(nonlinear_gaussian_data_function_SEXP)),
                load_matrix_function(nonlinear_gaussian_data_covariance_SEXP));
          } else {
            stop("nonlinear_gaussian_data specificiation is invalid.");
          }

          LikelihoodEstimator *new_likelihood_estimator =
              new ExactLikelihoodEstimator(rng_in, seed_in, data_in, proposal,
                                           factor_is_fixed_in_smc);

          full_index_vector.push_back(likelihood_estimators.size());
          without_cancelled_index_vector.push_back(
              likelihood_estimators.size());
          without_priors_index_vector.push_back(likelihood_estimators.size());

          if (factor_is_annealed) {
            likelihood_estimators.push_back(new AnnealedLikelihoodEstimator(
                rng_in, seed_in, data_in, new_likelihood_estimator, power,
                annealing_variable, false));
          } else {
            likelihood_estimators.push_back(new_likelihood_estimator);
          }
        } else if (current_factor.containsElementNamed(
                       "algorithmic_likelihood")) {
          if (Rf_isNewList(current_factor["algorithmic_likelihood"])) {
            List llhd_details = current_factor["algorithmic_likelihood"];

            if (!llhd_details.containsElementNamed("type")) {
              stop("algorithmic_likelihood must contain a type.");
            }

            if (!llhd_details.containsElementNamed("model")) {
              stop("algorithmic_likelihood must contain a model.");
            }

            if (!llhd_details.containsElementNamed("parameters")) {
              stop("algorithmic_likelihood must contain parameters.");
            }

            std::string llhd_type = llhd_details["type"];
            List llhd_model = llhd_details["model"];

            LikelihoodEstimator *new_likelihood_estimator;
            if (iequals(llhd_type, "kf")) {
              List kf_params =
                  get_kf_info(model_parameters, llhd_details, llhd_type);

              size_t iterations_to_store = kf_params[0];

              new_likelihood_estimator = get_kalman_filter(
                  data_in, llhd_model, model_parameters,
                  algorithm_parameter_list, iterations_to_store, false, "");
            } else if (iequals(llhd_type, "enkf")) {
              List enkf_params =
                  get_enkf_info(model_parameters, llhd_details, llhd_type);

              size_t number_of_ensemble_members = enkf_params[0];
              std::string shifter_name = enkf_params[1];
              size_t iterations_to_store = enkf_params[2];
              bool parallel = enkf_params[3];
              size_t grain_size = enkf_params[4];

              EnsembleShifter *shifter;
              if (iequals(shifter_name, "stochastic")) {
                shifter = new StochasticEnsembleShifter();
              } else if (iequals(shifter_name, "sqrt")) {
                // stop("Sqrt shifter not currently implemented for ensemble
                // Kalman filters.");
                //  To work, sqrt shifter needs construction of an "A" matrix (H
                //  matrix) that uses the state variables in the
                //  packing_instructions in EnsembleKalman and the state
                //  variables stored in the meascovest to construct this matrix.
                //  Can be done in the "setup" functions that are called before
                //  running the ensemble Kalman method.
                shifter = new SquareRootEnsembleShifter();
              } else if (iequals(shifter_name, "adjustment")) {
                shifter = new AdjustmentEnsembleShifter();
              } else {
                Rcout << "Shifter type: " << shifter_name << std::endl;
                Rcpp::stop("Invalid shifter type for EnKF");
              }

              new_likelihood_estimator = get_ensemble_kalman_filter(
                  rng_in, data_in, llhd_model, model_parameters,
                  algorithm_parameter_list, number_of_ensemble_members, shifter,
                  iterations_to_store, false, parallel, grain_size, "", seed_in,
                  data_created_in_get_measurement_covariance_estimators);
            } else if (iequals(llhd_type, "pf")) {
              List pf_params =
                  get_pf_info(model_parameters, llhd_details, llhd_type);

              size_t number_of_particles = pf_params[0];
              size_t iterations_to_store = pf_params[1];
              bool parallel = pf_params[2];
              size_t grain_size = pf_params[3];
              // an inner model that calls R cannot run in parallel threads
              if (parallel && llhd_model.containsElementNamed("uses_r_functions") &&
                  as<bool>(llhd_model["uses_r_functions"]))
                Rcpp::stop("A particle filter with parallel = TRUE needs a model written in C++ (or with ilike:: built-in functions): its model has functions written in R, Python or Julia, or R expressions, which cannot be called from parallel threads.");

              ParticleFilter *inner_pf = get_particle_filter(
                  rng_in, data_in, llhd_model, model_parameters,
                  algorithm_parameter_list, number_of_particles,
                  iterations_to_store, false, parallel, grain_size, "", seed_in,
                  data_created_in_get_likelihood_estimators,
                  data_created_in_get_measurement_covariance_estimators);
              // nested: keep its stored generations compactly and write them
              // concatenated over the outer points
              inner_pf->set_compact_stored_output(true);
              // threads for an inner particle filter with parallel = TRUE
              if (model.containsElementNamed("number_of_threads"))
                inner_pf->set_number_of_threads(size_t(as<double>(model["number_of_threads"])));
              new_likelihood_estimator = inner_pf;
            } else {
              stop("Invalid type for 'likelihood'.");
            }

            full_index_vector.push_back(likelihood_estimators.size());
            without_cancelled_index_vector.push_back(
                likelihood_estimators.size());
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
  if ((include_priors) && (any_annealing) &&
      (factors_affected_by_smc_sequence.empty())) {
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

std::vector<KalmanUpdater *> get_kalman_updaters(const List &model,
                                                 const List &model_parameters) {
  std::vector<KalmanUpdater *> updaters;

  if (model.containsElementNamed("factor")) {
    List factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];

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
        } else if (current_factor.containsElementNamed(
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

KalmanPredictor *get_kalman_predictor(const List &model,
                                      const List &model_parameters) {

  if (model.containsElementNamed("transition_model")) {
    List factors = model["transition_model"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];

        if (is_linear_gaussian_transition_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_transition_model(
              current_factor, model_parameters, false);
          if (!lg.is_function) {
            return new ExactKalmanPredictor(lg.matrix, lg.covariance);
          } else {
            return new ExactKalmanPredictor(lg.matrix_function,
                                            lg.covariance_function);
          }
        } else if (current_factor.containsElementNamed(
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

ProposalKernel *get_transition_proposal(const List &model,
                                        const List &model_parameters) {

  if (model.containsElementNamed("transition_proposal")) {
    List factors = model["transition_proposal"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];

        if (current_factor.containsElementNamed(
                "linear_gaussian_transition_covariance")) {
          if (!current_factor.containsElementNamed(
                  "linear_gaussian_transition_matrix")) {
            stop("linear_gaussian must contain both matrix and covariance.");
          }

          if (!current_factor.containsElementNamed("type")) {
            stop("linear_gaussian must contain a type.");
          }

          if (!current_factor.containsElementNamed(
                  "linear_gaussian_transition_variable")) {
            stop("linear_gaussian must contain "
                 "linear_gaussian_transition_variable.");
          }

          SEXP linear_gaussian_transition_matrix_SEXP =
              current_factor["linear_gaussian_transition_matrix"];
          SEXP linear_gaussian_transition_covariance_SEXP =
              current_factor["linear_gaussian_transition_covariance"];
          SEXP linear_gaussian_transition_variable_SEXP =
              current_factor["linear_gaussian_transition_variable"];

          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new LinearGaussianNoiseProposalKernel(
                load_string(linear_gaussian_transition_variable_SEXP),
                load_string(linear_gaussian_transition_variable_SEXP),
                load_matrix(linear_gaussian_transition_matrix_SEXP),
                load_matrix(linear_gaussian_transition_covariance_SEXP));
          } else if (type == 2) {
            proposal = new LinearGaussianNoiseFunctionProposalKernel(
                load_string(linear_gaussian_transition_variable_SEXP),
                load_string(linear_gaussian_transition_variable_SEXP),
                load_matrix_function(linear_gaussian_transition_matrix_SEXP),
                load_matrix_function(
                    linear_gaussian_transition_covariance_SEXP));
          } else {
            stop("linear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.containsElementNamed(
                       "nonlinear_gaussian_transition_covariance")) {

          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_transition_function")) {
            stop("nonlinear_gaussian_transition must contain both function and "
                 "covariance.");
          }

          if (!current_factor.containsElementNamed("type")) {
            stop("nonlinear_gaussian_transition must contain a type.");
          }

          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_transition_variable")) {
            stop("nonlinear_gaussian_transition must contain "
                 "nonlinear_gaussian_transition_variable.");
          }

          SEXP nonlinear_gaussian_transition_function_SEXP =
              current_factor["nonlinear_gaussian_transition_function"];
          SEXP nonlinear_gaussian_transition_covariance_SEXP =
              current_factor["nonlinear_gaussian_transition_covariance"];
          SEXP nonlinear_gaussian_transition_variable_SEXP =
              current_factor["nonlinear_gaussian_transition_variable"];
          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                load_string(nonlinear_gaussian_transition_variable_SEXP),
                std::make_shared<Transform>(load_transform(
                    nonlinear_gaussian_transition_function_SEXP)),
                load_matrix(nonlinear_gaussian_transition_covariance_SEXP));
          } else if (type == 2) {
            proposal = new NonLinearGaussianNoiseFunctionProposalKernel(
                load_string(nonlinear_gaussian_transition_variable_SEXP),
                std::make_shared<Transform>(load_transform(
                    nonlinear_gaussian_transition_function_SEXP)),
                load_matrix_function(
                    nonlinear_gaussian_transition_covariance_SEXP));
          } else {
            stop("nonlinear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.containsElementNamed("transition_proposal")) {

          if (Rf_isNewList(current_factor["transition_proposal"])) {

            List current_transition_proposal =
                current_factor["transition_proposal"];
            if (current_transition_proposal.containsElementNamed("type") &&
                current_transition_proposal.containsElementNamed("variables")) {
              std::string type = current_transition_proposal["type"];

              if (iequals(type, "norm_rw")) {
                List info = get_single_variable_one_double_parameter_info(
                    model_parameters, current_transition_proposal, type);

                std::string variable = info[0];
                double sd = info[1];

                return new GaussianRandomWalkProposalKernel(variable, sd);
              } else if (iequals(type, "mvnorm_rw")) {
                List info = get_single_variable_matrix_parameter_info(
                    model_parameters, current_transition_proposal, type);

                std::string variable = info[0];
                arma::mat cov = info[1];

                return new GaussianRandomWalkProposalKernel(variable, cov);
              }
              if (iequals(type, "unif_rw")) {
                List info = get_single_variable_one_double_parameter_info(
                    model_parameters, current_transition_proposal, type);

                std::string variable = info[0];
                double sd = info[1];

                return new UniformRandomWalkProposalKernel(variable, sd);
              }
              if (iequals(type, "mvunif_rw")) {
                List info = get_single_variable_vector_parameter_info(
                    model_parameters, current_transition_proposal, type);

                std::string variable = info[0];
                arma::colvec halfwidth = info[1];

                return new UniformRandomWalkProposalKernel(variable, halfwidth);
              } else {
                Rcout << "Transition model type " << type;
                stop("Transition model type unknown.");
              }
            } else {
              stop("Missing information in transition model section of model "
                   "file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.containsElementNamed(
                       "simulate_transition_proposal")) {
          SEXP simulate_transition_proposal_SEXP =
              current_factor["simulate_transition_proposal"];

          if (current_factor.containsElementNamed(
                  "evaluate_log_transition_proposal")) {
            SEXP evaluate_log_transition_proposal_SEXP =
                current_factor["evaluate_log_transition_proposal"];
            return new CustomGuidedNoParamsProposalKernel(
                load_simulate_guided_no_params_mcmc_proposal(
                    simulate_transition_proposal_SEXP),
                load_evaluate_log_guided_no_params_mcmc_proposal(
                    evaluate_log_transition_proposal_SEXP));
          } else {
            return new CustomGuidedNoParamsProposalKernel(
                load_simulate_guided_no_params_mcmc_proposal(
                    simulate_transition_proposal_SEXP));
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

ProposalKernel *get_transition_model(const List &model,
                                     const List &model_parameters) {

  if (model.containsElementNamed("transition_model")) {
    List factors = model["transition_model"];

    for (size_t i = 0; i < factors.size(); ++i) {

      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];

        if (is_linear_gaussian_transition_model(current_factor)) {
          LinearGaussianSpec lg = get_linear_gaussian_transition_model(
              current_factor, model_parameters, true);
          ProposalKernel *proposal;
          if (!lg.is_function) {
            proposal = new LinearGaussianNoiseProposalKernel(
                lg.variable, lg.variable, lg.matrix, lg.covariance);
          } else {
            proposal = new LinearGaussianNoiseFunctionProposalKernel(
                lg.variable, lg.variable, lg.matrix_function,
                lg.covariance_function);
          }

          return proposal;
        } else if (current_factor.containsElementNamed(
                       "nonlinear_gaussian_transition_covariance")) {

          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_transition_function")) {
            stop("nonlinear_gaussian_transition must contain both function and "
                 "covariance.");
          }

          if (!current_factor.containsElementNamed("type")) {
            stop("nonlinear_gaussian_transition must contain a type.");
          }

          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_transition_variable")) {
            stop("nonlinear_gaussian_transition must contain "
                 "nonlinear_gaussian_transition_variable.");
          }

          SEXP nonlinear_gaussian_transition_function_SEXP =
              current_factor["nonlinear_gaussian_transition_function"];
          SEXP nonlinear_gaussian_transition_covariance_SEXP =
              current_factor["nonlinear_gaussian_transition_covariance"];
          SEXP nonlinear_gaussian_transition_variable_SEXP =
              current_factor["nonlinear_gaussian_transition_variable"];
          size_t type = current_factor["type"];
          ProposalKernel *proposal;
          if (type == 1) {
            proposal = new NonLinearGaussianNoiseProposalKernel(
                load_string(nonlinear_gaussian_transition_variable_SEXP),
                std::make_shared<Transform>(load_transform(
                    nonlinear_gaussian_transition_function_SEXP)),
                load_matrix(nonlinear_gaussian_transition_covariance_SEXP));
          } else if (type == 2) {
            proposal = new NonLinearGaussianNoiseFunctionProposalKernel(
                load_string(nonlinear_gaussian_transition_variable_SEXP),
                std::make_shared<Transform>(load_transform(
                    nonlinear_gaussian_transition_function_SEXP)),
                load_matrix_function(
                    nonlinear_gaussian_transition_covariance_SEXP));
          } else {
            stop("nonlinear_gaussian_transition specificiation is invalid.");
          }

          return proposal;
        } else if (current_factor.containsElementNamed("transition_model")) {
          DistributionFactor *new_factor;

          if (Rf_isNewList(current_factor["transition_model"])) {

            List current_transition_model = current_factor["transition_model"];
            if (current_transition_model.containsElementNamed("type") &&
                current_transition_model.containsElementNamed("variables")) {
              std::string type = current_transition_model["type"];

              if (iequals(type, "norm_rw")) {
                List info = get_single_variable_one_double_parameter_info(
                    model_parameters, current_transition_model, type);

                std::string variable = info[0];
                double sd = info[1];

                return new GaussianRandomWalkProposalKernel(variable, sd);
              } else if (iequals(type, "mvnorm_rw")) {
                List info = get_single_variable_matrix_parameter_info(
                    model_parameters, current_transition_model, type);

                std::string variable = info[0];
                arma::mat cov = info[1];

                return new GaussianRandomWalkProposalKernel(variable, cov);
              }
              if (iequals(type, "unif_rw")) {
                List info = get_single_variable_one_double_parameter_info(
                    model_parameters, current_transition_model, type);

                std::string variable = info[0];
                double sd = info[1];

                return new UniformRandomWalkProposalKernel(variable, sd);
              }
              if (iequals(type, "mvunif_rw")) {
                List info = get_single_variable_vector_parameter_info(
                    model_parameters, current_transition_model, type);

                std::string variable = info[0];
                arma::colvec halfwidth = info[1];

                return new UniformRandomWalkProposalKernel(variable, halfwidth);
              } else {
                Rcout << "Transition model type " << type;
                stop("Transition model type unknown.");
              }
            } else {
              stop("Missing information in transition model section of model "
                   "file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.containsElementNamed(
                       "simulate_transition_model")) {
          SEXP simulate_transition_model_SEXP =
              current_factor["simulate_transition_model"];

          if (current_factor.containsElementNamed(
                  "evaluate_log_transition_model")) {
            SEXP evaluate_log_transition_model_SEXP =
                current_factor["evaluate_log_transition_model"];
            return new CustomNoParamsProposalKernel(
                load_simulate_no_params_mcmc_proposal(
                    simulate_transition_model_SEXP),
                load_evaluate_log_no_params_mcmc_proposal(
                    evaluate_log_transition_model_SEXP));
          } else {
            return new CustomNoParamsProposalKernel(
                load_simulate_no_params_mcmc_proposal(
                    simulate_transition_model_SEXP));
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

List get_prior_mean_and_covariance(const List &model,
                                   const List &model_parameters) {
  if (!model.containsElementNamed("factor"))
    stop(kf_prior_requirement + " No prior was found.");

  List factors = model["factor"];

  for (size_t i = 0; i < factors.size(); ++i) {

    if (!Rf_isNewList(factors[i]))
      stop("Error in factors section of model file.");

    List current_factor = factors[i];

    if (current_factor.containsElementNamed("evaluate_log_prior") ||
        current_factor.containsElementNamed("simulate_prior"))
      stop(kf_prior_requirement + " Factor " + std::to_string(i + 1) +
           " is a user-defined prior.");

    if (!current_factor.containsElementNamed("prior"))
      continue;

    if (!Rf_isNewList(current_factor["prior"]))
      stop("Error in factors section of model file.");

    List current_prior = current_factor["prior"];
    if (!(current_prior.containsElementNamed("type") &&
          current_prior.containsElementNamed("variables")))
      stop("Missing information for prior in model file.");

    std::string type = current_prior["type"];

    if (iequals(type, "norm")) {
      List info = get_single_variable_two_double_parameter_info(
          model_parameters, current_prior, type);

      double mean = info[1];
      double sd = info[2];

      arma::colvec mean_vec(1);
      mean_vec[0] = mean;

      arma::mat cov(1, 1);
      cov(0, 0) = sd * sd;

      return List::create(mean_vec, cov);
    } else if (iequals(type, "mvnorm")) {
      List info = get_single_variable_vector_and_matrix_parameter_info(
          model_parameters, current_prior, type);

      arma::colvec mean = info[1];
      arma::mat cov = info[2];

      return List::create(mean, cov);
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
    const List &model, const List &model_parameters,
    const std::vector<size_t> &factor_indices,
    std::shared_ptr<Transform> transform,
    std::vector<Data> &data_created_in_get_measurement_covariance_estimators) {

  std::vector<MeasurementCovarianceEstimator *>
      measurement_covariance_estimators;

  if (model.containsElementNamed("factor")) {
    List factors = model["factor"];

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

      if (Rf_isNewList(factors[indices_to_use[i]])) {
        List current_factor = factors[indices_to_use[i]];
        SimulateModelPtr simulate_data_model_in;
        if (current_factor.containsElementNamed("simulate_data_model")) {
          SEXP simulate_data_model_SEXP = current_factor["simulate_data_model"];
          simulate_data_model_in =
              load_simulate_data_model(simulate_data_model_SEXP);

          if (!current_factor.containsElementNamed("data_variable")) {
            stop("simulator-based model must be associated with a "
                 "data_variable.");
          }

          SEXP data_variable_SEXP = current_factor["data_variable"];

          std::string augmented_data_variable_names =
              load_string(data_variable_SEXP);
          std::vector<std::string> data_variable_names =
              split(augmented_data_variable_names, ';');

          MeasurementCovarianceEstimator *new_measurement_covariance_estimator;

          if (current_factor.containsElementNamed("summary_statistics")) {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                load_summary_statistics(summary_statistics_SEXP);

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

          if (current_factor.containsElementNamed("summary_statistics")) {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                load_summary_statistics(summary_statistics_SEXP);

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
        } else if (current_factor.containsElementNamed(
                       "nonlinear_gaussian_data_covariance")) {
          if (!current_factor.containsElementNamed(
                  "nonlinear_gaussian_data_function")) {
            stop("nonlinear_gaussian_data must contain both function and "
                 "covariance.");
          }

          if (!current_factor.containsElementNamed("type")) {
            stop("nonlinear_gaussian_data must contain a type.");
          }

          SEXP nonlinear_gaussian_data_function_SEXP =
              current_factor["nonlinear_gaussian_data_function"];
          SEXP nonlinear_gaussian_data_covariance_SEXP =
              current_factor["nonlinear_gaussian_data_covariance"];
          SEXP nonlinear_gaussian_data_variable_SEXP =
              current_factor["nonlinear_gaussian_data_variable"];
          size_t type = current_factor["type"];

          TransformPtr nonlinear_gaussian_data_function =
              load_transform(nonlinear_gaussian_data_function_SEXP);

          MeasurementCovarianceEstimator *new_measurement_covariance_estimator;

          if (current_factor.containsElementNamed("summary_statistics")) {
            SEXP summary_statistics_SEXP = current_factor["summary_statistics"];
            SummaryStatisticsPtr summary_stats =
                load_summary_statistics(summary_statistics_SEXP);

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
                      load_matrix(nonlinear_gaussian_data_covariance_SEXP),
                      load_string(nonlinear_gaussian_data_variable_SEXP));
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
                      load_matrix_function(
                          nonlinear_gaussian_data_covariance_SEXP),
                      load_string(nonlinear_gaussian_data_variable_SEXP));
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
                      load_matrix(nonlinear_gaussian_data_covariance_SEXP),
                      load_string(nonlinear_gaussian_data_variable_SEXP));
            } else if (type == 2) {
              new_measurement_covariance_estimator =
                  new DirectNonLinearGaussianMeasurementCovarianceEstimator(
                      rng_in, seed_in, data_in, transform, NULL,
                      std::make_shared<Transform>(
                          nonlinear_gaussian_data_function),
                      load_matrix_function(
                          nonlinear_gaussian_data_covariance_SEXP),
                      load_string(nonlinear_gaussian_data_variable_SEXP));
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
