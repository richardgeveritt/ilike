#ifndef TYPEDEFS
#define TYPEDEFS

// [[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"
#include "distributions.h"

// typedefs for function pointers

typedef double (*EvaluateLogDistributionPtr)(const Parameters &inputs);

typedef Parameters (*SimulateDistributionPtr)(RandomNumberGenerator &rng);

typedef double (*EvaluateLogLikelihoodPtr)(const Parameters &inputs, const Data &observed_data);

typedef double (*EstimateLogLikelihoodPtr)(const List &inputs, const List &observed_data, const List &auxiliary_variables);

typedef List (*SimulateModelPtr)(const List &inputs, const List &observed_data);

typedef List (*SimulateAuxiliaryVariablesPtr)(const List &inputs, const List &observed_data);

typedef XPtr<EstimateLogLikelihoodPtr> (*SetupLikelihoodEstimatorPtr)(const List &inputs, const List &auxiliary_variables);

typedef double (*EvaluateLogABCKernelPtr)(const arma::colvec &simulated_stats,
                const arma::colvec &observed_stats,
                const double &abc_tolerance);

typedef arma::colvec (*SummaryStatisticsPtr)(const List &observed_data);

typedef List (*GetDataFromSimulationPtr)(const List &simulation);

typedef Data (*DataPtr)(void);


EvaluateLogDistributionPtr load_evaluate_log_distribution(const SEXP &evaluate_log_distribution_SEXP);

SimulateDistributionPtr load_simulate_distribution(const SEXP &simulate_distribution_SEXP);

EvaluateLogLikelihoodPtr load_evaluate_log_likelihood(const SEXP &evaluate_log_likelihood_SEXP);

EstimateLogLikelihoodPtr load_estimate_log_likelihood(const SEXP &estimate_log_likelihood_SEXP);

SimulateModelPtr load_simulate_model(const SEXP &simulate_model_SEXP);

SimulateAuxiliaryVariablesPtr load_simulate_auxiliary_variables(const SEXP &simulate_auxiliary_variables_SEXP);

SetupLikelihoodEstimatorPtr load_setup_likelihood_estimator(const SEXP &setup_likelihood_estimator_SEXP);

EvaluateLogABCKernelPtr load_evaluate_log_abc_kernel(const SEXP &evaluate_log_abc_kernel_SEXP);

SummaryStatisticsPtr load_summary_statistics(const SEXP &summary_statistic_SEXP);

GetDataFromSimulationPtr load_get_data_from_simulation(const SEXP &get_data_from_simulation_SEXP);

//Data load_observed_data(const SEXP &observed_data_SEXP);

Data load_data(const SEXP &data_SEXP);

// List simulate_distribution_cpp(const SEXP &simulate_distribution_SEXP);
//
// double evaluate_log_distribution_cpp(const SEXP &evaluate_log_distribution_SEXP,
//                                      const List &parameter);
//
// double evaluate_log_likelihood_cpp(const SEXP &evaluate_log_likelihood_SEXP,
//                                    const List &parameter,
//                                    const List &observed_data);
//
// double estimate_log_likelihood_cpp(const SEXP &estimate_log_likelihood_SEXP,
//                                    const List &parameter,
//                                    const List &observed_data,
//                                    const List &auxiliary_variables);
//
// List simulate_model_cpp(const SEXP &simulate_model_SEXP,
//                         const List &parameter,
//                         const List &observed_data);
//
// List simulate_auxiliary_variables_cpp(const SEXP &simulate_auxiliary_variables_SEXP,
//                                       const List &parameter,
//                                       const List &observed_data);
//
// double evaluate_log_abc_kernel_cpp(const SEXP &evaluate_log_abc_kernel_SEXP,
//                                    const arma::colvec &simulated_data,
//                                    const arma::colvec &observed_data,
//                                    const double &abc_tolerance);
//
// arma::colvec summary_statistics_cpp(const SEXP &summary_statistics_SEXP,
//                                     const List &observed_data);
//
// List get_data_from_simulation_cpp(const SEXP &get_data_from_simulation_SEXP,
//                                            const List &simulation);

inline EvaluateLogDistributionPtr load_evaluate_log_distribution(const SEXP &evaluate_log_distribution_SEXP)
{
  XPtr<EvaluateLogDistributionPtr> evaluate_log_distribution_XPtr(evaluate_log_distribution_SEXP);
  return *evaluate_log_distribution_XPtr;
}

inline SimulateDistributionPtr load_simulate_distribution(const SEXP &simulate_distribution_SEXP)
{
  XPtr<SimulateDistributionPtr> simulate_distribution_XPtr(simulate_distribution_SEXP);
  return *simulate_distribution_XPtr;
}

inline EvaluateLogLikelihoodPtr load_evaluate_log_likelihood(const SEXP &evaluate_log_likelihood_SEXP)
{
  XPtr<EvaluateLogLikelihoodPtr> evaluate_log_likelihood_XPtr(evaluate_log_likelihood_SEXP);
  return *evaluate_log_likelihood_XPtr;
}

inline EstimateLogLikelihoodPtr load_estimate_log_likelihood(const SEXP &estimate_log_likelihood_SEXP)
{
  XPtr<EstimateLogLikelihoodPtr> estimate_log_likelihood_XPtr(estimate_log_likelihood_SEXP);
  return *estimate_log_likelihood_XPtr;
}

inline SimulateModelPtr load_simulate_model(const SEXP &simulate_model_SEXP)
{
  XPtr<SimulateModelPtr> simulate_model_XPtr(simulate_model_SEXP);
  return *simulate_model_XPtr;
}

inline SimulateAuxiliaryVariablesPtr load_simulate_auxiliary_variables(const SEXP &simulate_auxiliary_variables_SEXP)
{
  XPtr<SimulateAuxiliaryVariablesPtr> simulate_auxiliary_variables_XPtr(simulate_auxiliary_variables_SEXP);
  return *simulate_auxiliary_variables_XPtr;
}

inline SetupLikelihoodEstimatorPtr load_setup_likelihood_estimator(const SEXP &setup_likelihood_estimator_SEXP)
{
  XPtr<SetupLikelihoodEstimatorPtr> setup_likelihood_estimator_XPtr(setup_likelihood_estimator_SEXP);
  return *setup_likelihood_estimator_XPtr;
}

inline EvaluateLogABCKernelPtr load_evaluate_log_abc_kernel(const SEXP &evaluate_log_abc_kernel_SEXP)
{
  XPtr<EvaluateLogABCKernelPtr> evaluate_log_abc_kernel_XPtr(evaluate_log_abc_kernel_SEXP);
  return *evaluate_log_abc_kernel_XPtr;
}

inline SummaryStatisticsPtr load_summary_statistics(const SEXP &summary_statistics_SEXP)
{
  XPtr<SummaryStatisticsPtr> summary_statistics_XPtr(summary_statistics_SEXP);
  return *summary_statistics_XPtr;
}

inline GetDataFromSimulationPtr load_get_data_from_simulation(const SEXP &get_data_from_simulation_SEXP)
{
  XPtr<GetDataFromSimulationPtr> get_data_from_simulation_XPtr(get_data_from_simulation_SEXP);
  return *get_data_from_simulation_XPtr;
}

// inline Data load_observed_data(const SEXP &observed_data_SEXP)
// {
//   XPtr<Data> observed_data_XPtr(observed_data_SEXP);
//   return *observed_data_XPtr;
// }

inline Data load_data(const SEXP &data_SEXP)
{
  XPtr<DataPtr> data_XPtr(data_SEXP);
  return (*data_XPtr)();
}



// // [[Rcpp::export]]
// List simulate_distribution_cpp(const SEXP &simulate_distribution_SEXP)
// {
//   SimulateDistributionPtr simulate_distribution = load_simulate_distribution(simulate_distribution_SEXP);
//   return simulate_distribution();
// }
//
// // [[Rcpp::export]]
// double evaluate_log_distribution_cpp(const SEXP &evaluate_log_distribution_SEXP,
//                                      const List &parameter)
// {
//   EvaluateLogDistributionPtr evaluate_log_distribution = load_evaluate_log_distribution(evaluate_log_distribution_SEXP);
//   return evaluate_log_distribution(parameter);
// }
//
// // [[Rcpp::export]]
// double evaluate_log_likelihood_cpp(const SEXP &evaluate_log_likelihood_SEXP,
//                                    const List &parameter,
//                                    const List &observed_data)
// {
//   EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);
//   return evaluate_log_likelihood(parameter, observed_data);
// }
//
// // [[Rcpp::export]]
// double estimate_log_likelihood_cpp(const SEXP &estimate_log_likelihood_SEXP,
//                                    const List &parameter,
//                                    const List &observed_data,
//                                    const List &auxiliary_variables)
// {
//   EstimateLogLikelihoodPtr estimate_log_likelihood = load_estimate_log_likelihood(estimate_log_likelihood_SEXP);
//   return estimate_log_likelihood(parameter, observed_data, auxiliary_variables);
// }
//
// // [[Rcpp::export]]
// List simulate_model_cpp(const SEXP &simulate_model_SEXP,
//                         const List &parameter,
//                         const List &observed_data)
// {
//   SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);
//   return simulate_model(parameter, observed_data);
// }
//
// // [[Rcpp::export]]
// List simulate_auxiliary_variables_cpp(const SEXP &simulate_auxiliary_variables_SEXP,
//                                       const List &parameter,
//                                       const List &observed_data)
// {
//   SimulateModelPtr simulate_auxiliary_variables = load_simulate_auxiliary_variables(simulate_auxiliary_variables_SEXP);
//   return simulate_auxiliary_variables(parameter,observed_data);
// }
//
// // [[Rcpp::export]]
// double evaluate_log_abc_kernel_cpp(const SEXP &evaluate_log_abc_kernel_SEXP,
//                                    const arma::colvec &simulated_data,
//                                    const arma::colvec &observed_data,
//                                    const double &abc_tolerance)
// {
//   EvaluateLogABCKernelPtr evaluate_log_abc_kernel = load_evaluate_log_abc_kernel(evaluate_log_abc_kernel_SEXP);
//   return evaluate_log_abc_kernel(simulated_data, observed_data, abc_tolerance);
// }
//
// // [[Rcpp::export]]
// arma::colvec summary_statistics_cpp(const SEXP &summary_statistics_SEXP,
//                                     const List &observed_data)
// {
//   SummaryStatisticsPtr summary_statistics = load_summary_statistics(summary_statistics_SEXP);
//   return summary_statistics(observed_data);
// }
//
// // [[Rcpp::export]]
// List get_data_from_simulation_cpp(const SEXP &get_data_from_simulation_SEXP,
//                                            const List &simulation)
// {
//   GetDataFromSimulationPtr get_data_from_simulation = load_get_data_from_simulation(get_data_from_simulation_SEXP);
//   return get_data_from_simulation(simulation);
// }


#endif
