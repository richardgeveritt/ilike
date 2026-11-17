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
#include "custom_no_params_proposal_kernel.h"
#include "custom_proposal_kernel.h"
#include "adjustment_ensemble_shifter.h"
#include "always_smc_termination.h"
#include "annealed_likelihood_estimator.h"
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
void reset_informs();
bool ichar_equals(char a, char b);
bool iequals(const std::string &a, const std::string &b);
bool stob(const std::string &s);
std::vector<std::string> split(const std::string &str, const char delimiter);
bool isDoubleInList(const Rcpp::List &list, int index_from_one);
bool isStringInList(const Rcpp::List &list, int index_from_one);
bool isVectorInList(const Rcpp::List &list, int index_from_one);
bool isMatrixInList(const Rcpp::List &list, int index_from_one);
bool isNumericVectorInList(const Rcpp::List &list, int index_from_one);
bool isNumericMatrixInList(const Rcpp::List &list, int index_from_one);
bool isStringVectorInList(const Rcpp::List &list, const std::string &name);
bool isDoubleInList(const Rcpp::List &list, const std::string &name);
bool isNumericVectorInList(const Rcpp::List &list, const std::string &name);
bool isNumericMatrixInList(const Rcpp::List &list, const std::string &name);
Data get_data(const List &model);
bool extract_bool_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
int extract_int_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
double extract_double_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
std::string extract_string_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
arma::colvec extract_vector_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
arma::mat extract_matrix_parameter(const List &parameters_from_file, const List &model_parameters, size_t index);
List get_single_variable_one_double_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_two_double_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_vector_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_matrix_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_vector_and_matrix_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_matrix_and_double_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_single_variable_matrix_and_optional_double_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name, double default_double);
List get_single_variable_vector_and_matrix_and_double_parameter_info( const List &model_parameters, const List &current_distribution, const std::string &distribution_name);
List get_abc_euclidean_uniform_parameter_info(const List &model_parameters, const List &current_sbi, const std::string &sbi_name);
List get_abc_lp_uniform_parameter_info(const List &model_parameters, const List &current_sbi, const std::string &sbi_name);
List get_abc_enki_parameter_info(const List &model_parameters, const List &current_sbi, const std::string &sbi_name);
List get_sl_parameter_info(const List &model_parameters, const List &current_sbi, const std::string &sbi_name);
List get_kf_info(const List &model_parameters, const List &current_likelihood, const std::string &likelihood_name);
List get_enkf_info(const List &model_parameters, const List &current_likelihood, const std::string &likelihood_name);
List get_pf_info(const List &model_parameters, const List &current_likelihood, const std::string &likelihood_name);
std::vector<size_t> get_enk_likelihood_indices(const List &model, const List &model_parameters);
List get_ssm_info(const List &model, const List &model_parameters);
bool is_linear_gaussian_data_model(const List &current_factor);
bool is_linear_gaussian_transition_model(const List &current_factor);
List get_one_line_linear_gaussian(const List &current_factor, const std::string &block_name, size_t number_of_variables, const List &model_parameters);
LinearGaussianSpec get_linear_gaussian_data_model(const List &current_factor, const List &model_parameters);
LinearGaussianSpec get_linear_gaussian_transition_model(const List &current_factor, const List &model_parameters, bool variable_required);
std::string prior_cannot_be_annealed(size_t factor_index);
std::vector<LikelihoodEstimator *> get_likelihood_estimators( RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in, const List &model, const List &model_parameters, const List &algorithm_parameter_list, bool include_priors, const std::vector<std::string> &sequencer_types, const std::vector<std::string> &sequencer_variables, const std::vector<std::vector<double>> &sequencer_schedules, IndependentProposalKernel *proposal_in, VectorIndex *&without_cancelled_index, VectorIndex *&without_priors_index, VectorIndex *&full_index, bool &any_annealing, const std::vector<int> &factors_affected_by_smc_sequence, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
std::vector<KalmanUpdater *> get_kalman_updaters(const List &model, const List &model_parameters);
KalmanPredictor *get_kalman_predictor(const List &model, const List &model_parameters);
ProposalKernel *get_transition_proposal(const List &model, const List &model_parameters);
ProposalKernel *get_transition_model(const List &model, const List &model_parameters);
List get_prior_mean_and_covariance(const List &model, const List &model_parameters);
std::vector<MeasurementCovarianceEstimator *> get_measurement_covariance_estimators( RandomNumberGenerator *rng_in, size_t *seed_in, Data *data_in, const List &model, const List &model_parameters, const std::vector<size_t> &factor_indices, std::shared_ptr<Transform> transform, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
IndependentProposalKernel * get_proposal(const List &model, const List &model_parameters, Data *data);
IndependentProposalKernel *get_prior_as_proposal(const List &model, const List &model_parameters);
IndependentProposalKernel * get_prior_as_simulate_only_proposal(const List &model, const List &model_parameters);
ProposalKernel *get_mh_proposal(const List &current_proposal, const List &model_parameters, Data *data);
ProposalKernel *get_unadjusted_proposal(const List &current_proposal, const List &model_parameters, Data *data);
IndependentProposalKernel *get_imh_proposal(const List &current_proposal, const List &model_parameters, Data *data);
SymmetricProposalKernel *get_m_proposal(const List &current_proposal, const List &model_parameters, Data *data);
MCMCTermination *make_mcmc_termination(const List &model_parameters, const List &mcmc_termination_method);
void attach_mcmc_adaptor(ProposalKernel *proposal, const List &proposal_info, const std::string &prefix, const List &model_parameters);
void attach_smc_adaptor(ProposalKernel *proposal, const List &proposal_info, const List &model_parameters);
MCMC *make_mcmc(const List &model, const List &model_parameters, Data *data, const List &mcmc_termination_method, const List &mcmc_weights_method, Index *full_index);
std::vector<Parameters> make_initial_points(const List &initial_values);
size_t ilike_rdtsc();
Parameters make_algorithm_parameters(const List &algorithm_parameter_list);
SMCCriterion *get_resampling_method(const List &model, const List &model_parameters, size_t number_of_particles);
SMCCriterion *get_adaptive_target_method(const List &model, const List &model_parameters, const List &adaptive_target_method, size_t number_of_particles);
SMCTermination *get_smc_termination_method(const List &model, const List &model_parameters, const List &smc_termination_method);
EnsembleShifter *get_enk_shifter_method(const List &model, const List &model_parameters);
List get_smc_sequencer_info(const List &model, const List &model_parameters, const List &smc_sequencer_method);
ImportanceSampler *get_importance_sampler( RandomNumberGenerator *rng, Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_importance_points, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
double do_importance_sampler(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_importance_points, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
double do_importance_sampler_with_fixed_params( const List &model, const List &parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, size_t number_of_importance_points, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
void do_mcmc(const List &model, const List &parameters, const List &algorithm_parameter_list, const List &initial_values, const List &mcmc_termination_method, const List &mcmc_weights_method, size_t number_of_chains, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
void do_mcmc_with_fixed_params( const List &model, const List &parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, const List &initial_values, const List &mcmc_termination_method, const List &mcmc_weights_method, size_t number_of_chains, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
SMCMCMCMove *get_smc_mcmc_move( RandomNumberGenerator *rng, Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_particles, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &smc_sequencer_method, const List &adaptive_target_method, const List &smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
double do_smc_mcmc_move( const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_particles, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &smc_sequencer_method, const List &adaptive_target_method, const List &smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
double do_smc_mcmc_move_with_fixed_params( const List &model, const List &parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, size_t number_of_particles, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &smc_sequencer_method, const List &adaptive_target_method, const List &smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
EnsembleKalmanInversion *get_enki( RandomNumberGenerator *rng, Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_ensemble_members, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &enk_sequencer_method, const List &adaptive_target_method, const List &enk_termination_method, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
double do_enki(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_ensemble_members, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &enk_sequencer_method, const List &adaptive_target_method, const List &enk_termination_method, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
double do_enkmfds(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_particles, double Delta_t, const List &mcmc_termination_method, const List &mcmc_weights_method, const List &smc_sequencer_method, const List &adaptive_target_method, const List &smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
KalmanFilter *get_kalman_filter(Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t kf_iterations_to_store, bool write_to_file_at_each_iteration, const std::string &results_name_in);
double do_kalman_filter(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t kf_iterations_to_store, bool write_to_file_at_each_iteration, const String &results_name_in);
EnsembleKalmanFilter *get_ensemble_kalman_filter( RandomNumberGenerator *rng, Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_ensemble_members, EnsembleShifter *shifter, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
double do_ensemble_kalman_filter(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_ensemble_members, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
double do_ensemble_kalman_filter_with_fixed_params( const List &model, const List &parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, size_t number_of_ensemble_members, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
ParticleFilter *get_particle_filter( RandomNumberGenerator *rng, Data *the_data, const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_particles, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t *seed, std::vector<Data> &data_created_in_get_likelihood_estimators, std::vector<Data> &data_created_in_get_measurement_covariance_estimators);
double do_particle_filter(const List &model, const List &parameters, const List &algorithm_parameter_list, size_t number_of_particles, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
double do_particle_filter_with_fixed_params( const List &model, const List &parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, size_t number_of_particles, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String &results_name_in, size_t seed);
VectorIndex* user_factor_index(const IntegerVector &index_in, size_t number_of_factors);
double do_evaluate_log(const List &model, const List &model_parameters, const List &algorithm_parameter_list, const List &parameters_in, const IntegerVector &index_in, const String &results_name_in, size_t seed);
double do_evaluate_log_with_fixed_params(const List &model, const List &model_parameters, const List &algorithm_parameter_list, const List &fixed_parameter_list, const List &parameters_in, const IntegerVector &index_in, const String &results_name_in, size_t seed);
List shared_point_weights_check(const List &points_per_particle, const arma::colvec &normalised_log_weights);
List copy_audit(const List &model, const List &parameters, const List &algorithm_parameter_list, const List &mcmc_termination_method, const List &mcmc_weights_method, size_t seed_value);
void ilike_record_main_thread(DllInfo* dll);

#endif
