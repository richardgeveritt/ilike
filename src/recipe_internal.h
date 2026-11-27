#ifndef RECIPE_INTERNAL_H
#define RECIPE_INTERNAL_H

// Shared by the files that turn a recipe (the R list made by compile())
// into C++ objects and run algorithms on it: recipe_utilities.cpp (reading
// parameters and settings), recipe_factors.cpp (factors, likelihood
// estimators, transition models), recipe_mcmc.cpp (proposals, MCMC moves,
// SMC settings) and algorithm_interface.cpp (algorithms and the functions
// exported to R). dev/plan.md 3.5, CH4.

#include <set>
#include <algorithm> // std::equal
#include <cctype>    // std::tolower
#include <chrono>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "filesystem.h"
#include "ilike_hdf5_utils.h"
#include "algorithm_interface.h"
#include "builtin_registry.h"
#include "custom_no_params_proposal_kernel.h"
#include "custom_proposal_kernel.h"
#include "adjustment_ensemble_shifter.h"
#include "always_smc_termination.h"
#include "annealed_likelihood_estimator.h"
#include "cancellation_plan.h"
#include "unnormalised_normaliser_estimator.h"
#include "unnormalised_gradient_estimator.h"
#include "smc_path.h"
#include "barker_dynamics_proposal_kernel.h"
#include "cess_smc_criterion.h"
#include "composite_independent_proposal_kernel.h"
#include "custom_distribution_factor.h"
#include "custom_distribution_proposal_kernel.h"
#include "custom_guided_distribution_proposal_kernel.h"
#include "custom_guided_independent_proposal_kernel.h"
#include "custom_guided_no_params_proposal_kernel.h"
#include "custom_guided_no_params_symmetric_proposal_kernel.h"
#include "custom_guided_proposal_kernel.h"
#include "custom_guided_symmetric_proposal_kernel.h"
#include "custom_independent_proposal_kernel.h"
#include "custom_likelihood_factor.h"
#include "custom_no_params_symmetric_proposal_kernel.h"
#include "custom_symmetric_proposal_kernel.h"
#include "density_likelihood_estimator.h"
#include "deterministic_scan_mcmc.h"
#include "element_random_walk_proposal_kernel.h"
#include "in_place_element_scan_mcmc.h"
#include "element_adaptors.h"
#include "elementwise_transform_component.h"
#include "no_u_turn_mcmc.h"
#include "esjd_adaptor.h"
#include "transform.h"
#include "direct_gradient_estimator.h"
#include "direct_linear_gaussian_measurement_covariance_estimator.h"
#include "direct_nonlinear_gaussian_measurement_covariance_estimator.h"
#include "ensemble_kalman_filter.h"
#include "ensemble_kalman_inversion.h"
#include "ensemble_kalman_output.h"
#include "ess_smc_criterion.h"
#include "exact_kalman_predictor.h"
#include "exact_kalman_updater.h"
#include "exact_likelihood_estimator.h"
#include "gamma_distribution_factor.h"
#include "gamma_independent_proposal_kernel.h"
#include "gaussian_distribution_factor.h"
#include "gaussian_independent_proposal_kernel.h"
#include "gaussian_random_walk_proposal_kernel.h"
#include "gaussian_mcmc_adaptor.h"
#include "componentwise_random_walk_proposal_kernel.h"
#include "principal_components_proposal_kernel.h"
#include "proposal_parameter_adaptor.h"
#include "gaussian_smc_adaptor.h"
#include "generic_measurement_covariance_estimator.h"
#include "hmc_proposal_kernel.h"
#include "importance_sampler.h"
#include "iterations_mcmc_termination.h"
#include "kalman_filter.h"
#include "kalman_filter_output.h"
#include "langevin_proposal_kernel.h"
#include "likelihood_maker.h"
#include "chunk_parallel.h"
#include <memory>
#include <thread>
#include "linear_gaussian_noise_function_proposal_kernel.h"
#include "linear_gaussian_noise_proposal_kernel.h"
#include "loggaussian_distribution_factor.h"
#include "loggaussian_independent_proposal_kernel.h"
#include "measurement_covariance_estimator.h"
#include "metropolis_hastings_mcmc.h"
#include "metropolis_mcmc.h"
#include "mirror_proposal_kernel.h"
#include "mixed_generic_direct_gaussian_measurement_covariance_estimator.h"
#include "nonlinear_gaussian_noise_function_proposal_kernel.h"
#include "nonlinear_gaussian_noise_proposal_kernel.h"
#include "particle_filter.h"
#include "positive_smc_criterion.h"
#include "se_mcmc_termination.h"
#include "smc_mcmc_move.h"
#include "smc_output.h"
#include "square_root_ensemble_shifter.h"
#include "ensemble_kalman_mfds.h"
#include "stable_smc_termination.h"
#include "stochastic_ensemble_shifter.h"
#include "stochastic_scan_mcmc.h"
#include "transform.h"
#include "unadjusted_mcmc.h"
#include "uniform_random_walk_proposal_kernel.h"
#include "utils.h"
#include "shared_point_weights.h"
#include "vector_index.h"
#include "vector_factors.h"
#include "factor_variables.h"
#include "recipe_value.h"

using namespace Rcpp;
using namespace ilike;

// A linear-Gaussian model, given either by the separate blocks
// (linear_gaussian_{data,transition}_{matrix,covariance,variable,...}), with
// constant matrices (type 1) or matrix functions of the parameters (type 2),
// or by the one-line built-in
//   linear_gaussian_data_model,ilike::linear_gaussian_data_model(y,x,H,R)
//   linear_gaussian_transition_model,ilike::linear_gaussian_transition_model(x,A,Q)
// with constant matrices.
struct LinearGaussianSpec {
  std::string variable;
  std::string state_variable_string;
  std::vector<std::string> state_variables;
  bool is_function = false;
  arma::mat matrix;
  arma::mat covariance;
  GetMatrixPtr matrix_function = NULL;
  GetMatrixPtr covariance_function = NULL;
};
// The Kalman filter needs the prior on the initial state to be Gaussian with
// a known mean and covariance, so it must be a built-in Gaussian prior.
const std::string kf_prior_requirement =
    "The Kalman filter requires the prior on the initial state to be a "
    "built-in Gaussian: prior,ilike::norm(variable,mean,sd) for a "
    "one-dimensional state or prior,ilike::mvnorm(variable,mean,covariance).";
// Reads the optional "<prefix>_factor_index" (factors used in the move's
// acceptance ratio) and "<prefix>_cut_factor_index" (factors that depend on
// what the move changes but are deliberately left out, and are re-evaluated
// after acceptance) from one proposal block of the recipe, and sets them on
// the move. Indices are 1-based in the recipe. Defaults: no factor index means
// all factors except the cut ones; no cut index means no cut. Factors in
// neither set are assumed to cancel in the acceptance ratio.
// The factors one move uses and cuts (0-based), with defaults filled in, and a
// label for error messages.
struct MoveFactors {
  std::string label;
  bool unadjusted;
  std::vector<size_t> used;
  std::vector<size_t> cut;
};

void inform_once(const std::string &note);
bool ichar_equals(char a, char b);
bool iequals(const std::string &a, const std::string &b);
std::vector<std::string> split(const std::string &str, const char delimiter);
bool isDoubleInList(const RecipeValue &list, const std::string &name);
bool isNumericVectorInList(const RecipeValue &list, const std::string &name);
bool isNumericMatrixInList(const RecipeValue &list, const std::string &name);
Data get_data(const RecipeValue &model);
RecipeValue get_abc_euclidean_uniform_parameter_info(const RecipeValue &model_parameters, const RecipeValue &current_sbi, const std::string &sbi_name);
RecipeValue get_abc_lp_uniform_parameter_info(const RecipeValue &model_parameters, const RecipeValue &current_sbi, const std::string &sbi_name);
RecipeValue get_abc_enki_parameter_info(const RecipeValue &model_parameters, const RecipeValue &current_sbi, const std::string &sbi_name);
RecipeValue get_sl_parameter_info(const RecipeValue &model_parameters, const RecipeValue &current_sbi, const std::string &sbi_name);
RecipeValue get_kf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood, const std::string &likelihood_name);
RecipeValue get_enkf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood, const std::string &likelihood_name);
RecipeValue get_pf_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood, const std::string &likelihood_name);
RecipeValue get_is_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood);
RecipeValue get_integrate_smc_info(const RecipeValue &model_parameters, const RecipeValue &current_likelihood);
std::vector<size_t> get_enk_likelihood_indices(const RecipeValue &model, const RecipeValue &model_parameters);
RecipeValue get_ssm_info(const RecipeValue &model, const RecipeValue &model_parameters);
bool is_linear_gaussian_data_model(const RecipeValue &current_factor);
bool is_linear_gaussian_transition_model(const RecipeValue &current_factor);
RecipeValue get_one_line_linear_gaussian(const RecipeValue &current_factor, const std::string &block_name, size_t number_of_variables, const RecipeValue &model_parameters);
LinearGaussianSpec get_linear_gaussian_data_model(const RecipeValue &current_factor, const RecipeValue &model_parameters);
LinearGaussianSpec get_linear_gaussian_transition_model(const RecipeValue &current_factor, const RecipeValue &model_parameters, bool variable_required);
TargetSpecification target_specification( const RecipeValue &model, const std::vector<std::string> &sequencer_types, const std::vector<std::string> &sequencer_variables, const std::vector<int> &factors_affected_by_smc_sequence);
std::vector<LikelihoodEstimator *> get_likelihood_estimators( RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in, const RecipeValue &model, const RecipeValue &model_parameters, const Parameters &algorithm_parameters, const KernelDeclaration &importance_proposal, const std::vector<std::string> &sequencer_types, const std::vector<std::string> &sequencer_variables, const std::vector<std::vector<double>> &sequencer_schedules, IndependentProposalKernel *proposal_in, VectorIndex *&without_cancelled_index, VectorIndex *&without_priors_index, VectorIndex *&full_index, CancellationPlan &plan, const std::vector<int> &factors_affected_by_smc_sequence, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators, const DataTemperingOptions &data_tempering = DataTemperingOptions());
std::vector<KalmanUpdater *> get_kalman_updaters(const RecipeValue &model, const RecipeValue &model_parameters);
KalmanPredictor *get_kalman_predictor(const RecipeValue &model, const RecipeValue &model_parameters);
ProposalKernel *get_transition_proposal(const RecipeValue &model, const RecipeValue &model_parameters);
ProposalKernel *get_transition_model(const RecipeValue &model, const RecipeValue &model_parameters);
RecipeValue get_prior_mean_and_covariance(const RecipeValue &model, const RecipeValue &model_parameters);
// The means and covariances of a prior made only of built-in Gaussian or
// log-normal priors (norm, mvnorm, lnorm, mvlnorm), one entry per variable, on
// the Gaussian scale; log_variables are the log-normal ones (MFDS, dev/plan.md 3.16).
void get_gaussian_prior_moments(const RecipeValue &model, const RecipeValue &model_parameters,
                                Parameters &means, Parameters &covariances,
                                std::vector<std::string> &log_variables);
std::vector<MeasurementCovarianceEstimator *> get_measurement_covariance_estimators( RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in, const RecipeValue &model, const RecipeValue &model_parameters, const std::vector<size_t> &factor_indices, std::shared_ptr<Transform> transform, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
IndependentProposalKernel * get_proposal(const RecipeValue &model, const RecipeValue &model_parameters, Data *data);
IndependentProposalKernel *get_prior_as_proposal(const RecipeValue &model, const RecipeValue &model_parameters);
IndependentProposalKernel * get_prior_as_simulate_only_proposal(const RecipeValue &model, const RecipeValue &model_parameters);
ProposalKernel *get_mh_proposal(const RecipeValue &current_proposal, const RecipeValue &model_parameters, Data *data);
ProposalKernel *get_unadjusted_proposal(const RecipeValue &current_proposal, const RecipeValue &model_parameters, Data *data);
IndependentProposalKernel *get_imh_proposal(const RecipeValue &current_proposal, const RecipeValue &model_parameters, Data *data);
SymmetricProposalKernel *get_m_proposal(const RecipeValue &current_proposal, const RecipeValue &model_parameters, Data *data);
MCMCTermination *make_mcmc_termination(const RecipeValue &model_parameters, const RecipeValue &mcmc_termination_method);
void attach_mcmc_adaptor(ProposalKernel *proposal, const RecipeValue &proposal_info, const std::string &prefix, const RecipeValue &model_parameters);
// the built-in transforms from the recipe (C209); NULL if none
Transform *transform_from_recipe(const RecipeValue &list);
void attach_smc_adaptor(ProposalKernel *proposal, const RecipeValue &proposal_info, const RecipeValue &model_parameters);
MCMC *make_mcmc(const RecipeValue &model, const RecipeValue &model_parameters, Data *data, const RecipeValue &mcmc_termination_method, const RecipeValue &mcmc_weights_method, Index *full_index);
SMCCriterion *get_resampling_method(const RecipeValue &model, const RecipeValue &model_parameters, size_t number_of_particles);
SMCCriterion *get_adaptive_target_method(const RecipeValue &model, const RecipeValue &model_parameters, const RecipeValue &adaptive_target_method, size_t number_of_particles);
SMCTermination *get_smc_termination_method(const RecipeValue &model, const RecipeValue &model_parameters, const RecipeValue &smc_termination_method);
EnsembleShifter *get_enk_shifter_method(const RecipeValue &model, const RecipeValue &model_parameters);
// An importance sampler or SMC inside a factor (integrate u with is(…) /
// smc(…), dev/plan.md 3.21 (b), (c2)): the outer SMC sequences that change the factor (none if it
// is fixed along the outer sequence). The inner factors then change with the
// outer target, and the inner points are kept and weighted again under each
// new target (Kerama, Thorne & Everitt's ABC-SMC, arXiv:2211.02172, alg. 3).
struct NestedSequence {
  std::vector<std::string> types;
  std::vector<std::string> variables;
};
RecipeValue get_smc_sequencer_info(const RecipeValue &model, const RecipeValue &model_parameters, const RecipeValue &smc_sequencer_method,
                                   BuiltinRole role = BuiltinRole::SMCSequence);
ImportanceSampler *get_importance_sampler( RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t number_of_importance_points, bool parallel_in, size_t grain_size_in, const std::string &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators, const NestedSequence *nested = NULL);
SMCMCMCMove *get_smc_mcmc_move( RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t number_of_particles, const RecipeValue &mcmc_termination_method, const RecipeValue &mcmc_weights_method, const RecipeValue &smc_sequencer_method, const RecipeValue &adaptive_target_method, const RecipeValue &smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const std::string &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators, const NestedSequence *nested = NULL);
EnsembleKalmanInversion *get_enki( RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t number_of_ensemble_members, const RecipeValue &mcmc_termination_method, const RecipeValue &mcmc_weights_method, const RecipeValue &enk_sequencer_method, const RecipeValue &adaptive_target_method, const RecipeValue &enk_termination_method, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const std::string &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
KalmanFilter *get_kalman_filter(Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t kf_iterations_to_store, bool write_to_file_at_each_iteration, const std::string &results_name_in);
EnsembleKalmanFilter *get_ensemble_kalman_filter( RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t number_of_ensemble_members, EnsembleShifter *shifter, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const std::string &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
ParticleFilter *get_particle_filter( RandomNumberGenerator *rng, Data *the_data, const RecipeValue &model, const RecipeValue &parameters, const Parameters &algorithm_parameters, size_t number_of_particles, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const std::string &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
void ilike_record_main_thread(DllInfo* dll);
// (dev/plan.md 3.22 U4) the general engine of an unnormalised distribution written as a block
SMCPath *make_smc_path(RandomNumberGenerator *rng, size_t *seed, Data *data, const RecipeValue &model, const RecipeValue &model_parameters, const std::string &x, const std::string &mode, size_t particles, size_t steps, size_t sweeps, double resampling, size_t waste_free);

#endif
