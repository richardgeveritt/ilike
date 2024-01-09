// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/ilike.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ilike_rdtsc
size_t ilike_rdtsc();
RcppExport SEXP _ilike_ilike_rdtsc() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(ilike_rdtsc());
    return rcpp_result_gen;
END_RCPP
}
// do_importance_sampler
double do_importance_sampler(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_importance_points, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_importance_sampler(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_importance_pointsSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_importance_points(number_of_importance_pointsSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_importance_sampler(model, parameters, algorithm_parameter_list, number_of_importance_points, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}
// do_mcmc
void do_mcmc(const List& model, const List& parameters, const List& algorithm_parameter_list, const List& initial_values, const List& mcmc_termination_method, const List& mcmc_weights_method, size_t number_of_chains, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_mcmc(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP initial_valuesSEXP, SEXP mcmc_termination_methodSEXP, SEXP mcmc_weights_methodSEXP, SEXP number_of_chainsSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< const List& >::type initial_values(initial_valuesSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_termination_method(mcmc_termination_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_weights_method(mcmc_weights_methodSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_chains(number_of_chainsSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    do_mcmc(model, parameters, algorithm_parameter_list, initial_values, mcmc_termination_method, mcmc_weights_method, number_of_chains, parallel_in, grain_size_in, results_name_in, seed);
    return R_NilValue;
END_RCPP
}
// do_smc_mcmc_move
double do_smc_mcmc_move(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_particles, const List& mcmc_termination_method, const List& mcmc_weights_method, const List& smc_sequencer_method, const List& adaptive_target_method, const List& smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_smc_mcmc_move(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_particlesSEXP, SEXP mcmc_termination_methodSEXP, SEXP mcmc_weights_methodSEXP, SEXP smc_sequencer_methodSEXP, SEXP adaptive_target_methodSEXP, SEXP smc_termination_methodSEXP, SEXP smc_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_particles(number_of_particlesSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_termination_method(mcmc_termination_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_weights_method(mcmc_weights_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type smc_sequencer_method(smc_sequencer_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type adaptive_target_method(adaptive_target_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type smc_termination_method(smc_termination_methodSEXP);
    Rcpp::traits::input_parameter< size_t >::type smc_iterations_to_store(smc_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_smc_mcmc_move(model, parameters, algorithm_parameter_list, number_of_particles, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}
// do_enki
double do_enki(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_ensemble_members, const List& mcmc_termination_method, const List& mcmc_weights_method, const List& enk_sequencer_method, const List& adaptive_target_method, const List& enk_termination_method, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_enki(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_ensemble_membersSEXP, SEXP mcmc_termination_methodSEXP, SEXP mcmc_weights_methodSEXP, SEXP enk_sequencer_methodSEXP, SEXP adaptive_target_methodSEXP, SEXP enk_termination_methodSEXP, SEXP enk_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_ensemble_members(number_of_ensemble_membersSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_termination_method(mcmc_termination_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_weights_method(mcmc_weights_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type enk_sequencer_method(enk_sequencer_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type adaptive_target_method(adaptive_target_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type enk_termination_method(enk_termination_methodSEXP);
    Rcpp::traits::input_parameter< size_t >::type enk_iterations_to_store(enk_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_enki(model, parameters, algorithm_parameter_list, number_of_ensemble_members, mcmc_termination_method, mcmc_weights_method, enk_sequencer_method, adaptive_target_method, enk_termination_method, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}
// do_enkmfds
double do_enkmfds(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_particles, double Delta_t, const List& mcmc_termination_method, const List& mcmc_weights_method, const List& smc_sequencer_method, const List& adaptive_target_method, const List& smc_termination_method, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_enkmfds(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_particlesSEXP, SEXP Delta_tSEXP, SEXP mcmc_termination_methodSEXP, SEXP mcmc_weights_methodSEXP, SEXP smc_sequencer_methodSEXP, SEXP adaptive_target_methodSEXP, SEXP smc_termination_methodSEXP, SEXP smc_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_particles(number_of_particlesSEXP);
    Rcpp::traits::input_parameter< double >::type Delta_t(Delta_tSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_termination_method(mcmc_termination_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type mcmc_weights_method(mcmc_weights_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type smc_sequencer_method(smc_sequencer_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type adaptive_target_method(adaptive_target_methodSEXP);
    Rcpp::traits::input_parameter< const List& >::type smc_termination_method(smc_termination_methodSEXP);
    Rcpp::traits::input_parameter< size_t >::type smc_iterations_to_store(smc_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_enkmfds(model, parameters, algorithm_parameter_list, number_of_particles, Delta_t, mcmc_termination_method, mcmc_weights_method, smc_sequencer_method, adaptive_target_method, smc_termination_method, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}
// do_kalman_filter
double do_kalman_filter(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t kf_iterations_to_store, bool write_to_file_at_each_iteration, const String& results_name_in);
RcppExport SEXP _ilike_do_kalman_filter(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP kf_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP results_name_inSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type kf_iterations_to_store(kf_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    rcpp_result_gen = Rcpp::wrap(do_kalman_filter(model, parameters, algorithm_parameter_list, kf_iterations_to_store, write_to_file_at_each_iteration, results_name_in));
    return rcpp_result_gen;
END_RCPP
}
// do_ensemble_kalman_filter
double do_ensemble_kalman_filter(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_ensemble_members, size_t enk_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_ensemble_kalman_filter(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_ensemble_membersSEXP, SEXP enk_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_ensemble_members(number_of_ensemble_membersSEXP);
    Rcpp::traits::input_parameter< size_t >::type enk_iterations_to_store(enk_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_ensemble_kalman_filter(model, parameters, algorithm_parameter_list, number_of_ensemble_members, enk_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}
// do_particle_filter
double do_particle_filter(const List& model, const List& parameters, const List& algorithm_parameter_list, size_t number_of_particles, size_t smc_iterations_to_store, bool write_to_file_at_each_iteration, bool parallel_in, size_t grain_size_in, const String& results_name_in, size_t seed);
RcppExport SEXP _ilike_do_particle_filter(SEXP modelSEXP, SEXP parametersSEXP, SEXP algorithm_parameter_listSEXP, SEXP number_of_particlesSEXP, SEXP smc_iterations_to_storeSEXP, SEXP write_to_file_at_each_iterationSEXP, SEXP parallel_inSEXP, SEXP grain_size_inSEXP, SEXP results_name_inSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type model(modelSEXP);
    Rcpp::traits::input_parameter< const List& >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< const List& >::type algorithm_parameter_list(algorithm_parameter_listSEXP);
    Rcpp::traits::input_parameter< size_t >::type number_of_particles(number_of_particlesSEXP);
    Rcpp::traits::input_parameter< size_t >::type smc_iterations_to_store(smc_iterations_to_storeSEXP);
    Rcpp::traits::input_parameter< bool >::type write_to_file_at_each_iteration(write_to_file_at_each_iterationSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel_in(parallel_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type grain_size_in(grain_size_inSEXP);
    Rcpp::traits::input_parameter< const String& >::type results_name_in(results_name_inSEXP);
    Rcpp::traits::input_parameter< size_t >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(do_particle_filter(model, parameters, algorithm_parameter_list, number_of_particles, smc_iterations_to_store, write_to_file_at_each_iteration, parallel_in, grain_size_in, results_name_in, seed));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ilike_ilike_rdtsc", (DL_FUNC) &_ilike_ilike_rdtsc, 0},
    {"_ilike_do_importance_sampler", (DL_FUNC) &_ilike_do_importance_sampler, 8},
    {"_ilike_do_mcmc", (DL_FUNC) &_ilike_do_mcmc, 11},
    {"_ilike_do_smc_mcmc_move", (DL_FUNC) &_ilike_do_smc_mcmc_move, 15},
    {"_ilike_do_enki", (DL_FUNC) &_ilike_do_enki, 15},
    {"_ilike_do_enkmfds", (DL_FUNC) &_ilike_do_enkmfds, 16},
    {"_ilike_do_kalman_filter", (DL_FUNC) &_ilike_do_kalman_filter, 6},
    {"_ilike_do_ensemble_kalman_filter", (DL_FUNC) &_ilike_do_ensemble_kalman_filter, 10},
    {"_ilike_do_particle_filter", (DL_FUNC) &_ilike_do_particle_filter, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_ilike(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
