# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ilike_rdtsc <- function() {
    .Call(`_ilike_ilike_rdtsc`)
}

do_importance_sampler <- function(model, parameters, algorithm_parameter_list, number_of_importance_points, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_importance_sampler`, model, parameters, algorithm_parameter_list, number_of_importance_points, parallel_in, grain_size_in, results_name_in, seed)
}

do_importance_sampler_with_fixed_params <- function(model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_importance_points, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_importance_sampler_with_fixed_params`, model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_importance_points, parallel_in, grain_size_in, results_name_in, seed)
}

do_mcmc <- function(model, parameters, algorithm_parameter_list, initial_values, mcmc_termination_method, mcmc_weights_method, number_of_chains, parallel_in, grain_size_in, results_name_in, seed) {
    invisible(.Call(`_ilike_do_mcmc`, model, parameters, algorithm_parameter_list, initial_values, mcmc_termination_method, mcmc_weights_method, number_of_chains, parallel_in, grain_size_in, results_name_in, seed))
}

do_mcmc_with_fixed_params <- function(model, parameters, algorithm_parameter_list, fixed_parameter_list, initial_values, mcmc_termination_method, mcmc_weights_method, number_of_chains, parallel_in, grain_size_in, results_name_in, seed) {
    invisible(.Call(`_ilike_do_mcmc_with_fixed_params`, model, parameters, algorithm_parameter_list, fixed_parameter_list, initial_values, mcmc_termination_method, mcmc_weights_method, number_of_chains, parallel_in, grain_size_in, results_name_in, seed))
}

do_smc_mcmc_move <- function(model, parameters, algorithm_parameter_list, number_of_particles, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_smc_mcmc_move`, model, parameters, algorithm_parameter_list, number_of_particles, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_enki <- function(model, parameters, algorithm_parameter_list, number_of_ensemble_members, mcmc_termination_method, mcmc_weights_method, enk_sequencer_method, adaptive_target_method, enk_termination_method, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_enki`, model, parameters, algorithm_parameter_list, number_of_ensemble_members, mcmc_termination_method, mcmc_weights_method, enk_sequencer_method, adaptive_target_method, enk_termination_method, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_enkmfds <- function(model, parameters, algorithm_parameter_list, number_of_particles, Delta_t, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_enkmfds`, model, parameters, algorithm_parameter_list, number_of_particles, Delta_t, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_kalman_filter <- function(model, parameters, algorithm_parameter_list, kf_iterations_to_store, write_to_file_at_each_iteration, results_name_in) {
    .Call(`_ilike_do_kalman_filter`, model, parameters, algorithm_parameter_list, kf_iterations_to_store, write_to_file_at_each_iteration, results_name_in)
}

do_ensemble_kalman_filter <- function(model, parameters, algorithm_parameter_list, number_of_ensemble_members, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_ensemble_kalman_filter`, model, parameters, algorithm_parameter_list, number_of_ensemble_members, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_ensemble_kalman_filter_with_fixed_params <- function(model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_ensemble_members, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_ensemble_kalman_filter_with_fixed_params`, model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_ensemble_members, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_particle_filter <- function(model, parameters, algorithm_parameter_list, number_of_particles, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_particle_filter`, model, parameters, algorithm_parameter_list, number_of_particles, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

do_particle_filter_with_fixed_params <- function(model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_particles, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed) {
    .Call(`_ilike_do_particle_filter_with_fixed_params`, model, parameters, algorithm_parameter_list, fixed_parameter_list, number_of_particles, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed)
}

