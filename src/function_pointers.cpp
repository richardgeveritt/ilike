#include "function_pointers.h"

EvaluateLogDistributionPtr load_evaluate_log_distribution(const SEXP &evaluate_log_distribution_SEXP)
{
  XPtr<EvaluateLogDistributionPtr> evaluate_log_distribution_XPtr(evaluate_log_distribution_SEXP);
  return *evaluate_log_distribution_XPtr;
}

SimulateDistributionPtr load_simulate_distribution(const SEXP &simulate_distribution_SEXP)
{
  XPtr<SimulateDistributionPtr> simulate_distribution_XPtr(simulate_distribution_SEXP);
  return *simulate_distribution_XPtr;
}

EvaluateLogLikelihoodPtr load_evaluate_log_likelihood(const SEXP &evaluate_log_likelihood_SEXP)
{
  XPtr<EvaluateLogLikelihoodPtr> evaluate_log_likelihood_XPtr(evaluate_log_likelihood_SEXP);
  return *evaluate_log_likelihood_XPtr;
}

EstimateLogLikelihoodPtr load_estimate_log_likelihood(const SEXP &estimate_log_likelihood_SEXP)
{
  XPtr<EstimateLogLikelihoodPtr> estimate_log_likelihood_XPtr(estimate_log_likelihood_SEXP);
  return *estimate_log_likelihood_XPtr;
}

SimulateModelPtr load_simulate_model(const SEXP &simulate_model_SEXP)
{
  XPtr<SimulateModelPtr> simulate_model_XPtr(simulate_model_SEXP);
  return *simulate_model_XPtr;
}

SimulateAuxiliaryVariablesPtr load_simulate_auxiliary_variables(const SEXP &simulate_auxiliary_variables_SEXP)
{
  XPtr<SimulateAuxiliaryVariablesPtr> simulate_auxiliary_variables_XPtr(simulate_auxiliary_variables_SEXP);
  return *simulate_auxiliary_variables_XPtr;
}

SetupLikelihoodEstimatorPtr load_setup_likelihood_estimator(const SEXP &setup_likelihood_estimator_SEXP)
{
  XPtr<SetupLikelihoodEstimatorPtr> setup_likelihood_estimator_XPtr(setup_likelihood_estimator_SEXP);
  return *setup_likelihood_estimator_XPtr;
}

EvaluateLogABCKernelPtr load_evaluate_log_abc_kernel(const SEXP &evaluate_log_abc_kernel_SEXP)
{
  XPtr<EvaluateLogABCKernelPtr> evaluate_log_abc_kernel_XPtr(evaluate_log_abc_kernel_SEXP);
  return *evaluate_log_abc_kernel_XPtr;
}

SummaryStatisticsPtr load_summary_statistics(const SEXP &summary_statistics_SEXP)
{
  XPtr<SummaryStatisticsPtr> summary_statistics_XPtr(summary_statistics_SEXP);
  return *summary_statistics_XPtr;
}

GetDataFromSimulationPtr load_get_data_from_simulation(const SEXP &get_data_from_simulation_SEXP)
{
  XPtr<GetDataFromSimulationPtr> get_data_from_simulation_XPtr(get_data_from_simulation_SEXP);
  return *get_data_from_simulation_XPtr;
}



// [[Rcpp::export]]
NumericVector simulate_distribution_cpp(const SEXP &simulate_distribution_SEXP)
{
  SimulateDistributionPtr simulate_distribution = load_simulate_distribution(simulate_distribution_SEXP);
  return simulate_distribution();
}

// [[Rcpp::export]]
double evaluate_log_distribution_cpp(const SEXP &evaluate_log_distribution_SEXP,
                                     const NumericVector &parameter)
{
  EvaluateLogDistributionPtr evaluate_log_distribution = load_evaluate_log_distribution(evaluate_log_distribution_SEXP);
  return evaluate_log_distribution(parameter);
}

// [[Rcpp::export]]
double evaluate_log_likelihood_cpp(const SEXP &evaluate_log_likelihood_SEXP,
                                   const NumericVector &parameter,
                                   const NumericMatrix &data)
{
  EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);
  return evaluate_log_likelihood(parameter, data);
}

// [[Rcpp::export]]
double estimate_log_likelihood_cpp(const SEXP &estimate_log_likelihood_SEXP,
                                   const NumericVector &parameter,
                                   const NumericMatrix &data,
                                   const List &auxiliary_variables)
{
  EstimateLogLikelihoodPtr estimate_log_likelihood = load_estimate_log_likelihood(estimate_log_likelihood_SEXP);
  return estimate_log_likelihood(parameter, data, auxiliary_variables);
}

// [[Rcpp::export]]
List simulate_model_cpp(const SEXP &simulate_model_SEXP,
                        const NumericVector &parameter,
                        const NumericMatrix &data)
{
  SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);
  return simulate_model(parameter, data);
}

// [[Rcpp::export]]
List simulate_auxiliary_variables_cpp(const SEXP &simulate_auxiliary_variables_SEXP,
                                      const NumericVector &parameter,
                                      const NumericMatrix &data)
{
  SimulateModelPtr simulate_auxiliary_variables = load_simulate_auxiliary_variables(simulate_auxiliary_variables_SEXP);
  return simulate_auxiliary_variables(parameter,data);
}

// [[Rcpp::export]]
double evaluate_log_abc_kernel_cpp(const SEXP &evaluate_log_abc_kernel_SEXP,
                                   const arma::colvec &simulated_data,
                                   const arma::colvec &data,
                                   const double &abc_tolerance)
{
  EvaluateLogABCKernelPtr evaluate_log_abc_kernel = load_evaluate_log_abc_kernel(evaluate_log_abc_kernel_SEXP);
  return evaluate_log_abc_kernel(simulated_data, data, abc_tolerance);
}

// [[Rcpp::export]]
arma::colvec summary_statistics_cpp(const SEXP &summary_statistics_SEXP,
                                    const NumericMatrix &data)
{
  SummaryStatisticsPtr summary_statistics = load_summary_statistics(summary_statistics_SEXP);
  return summary_statistics(data);
}

// [[Rcpp::export]]
NumericMatrix get_data_from_simulation_cpp(const SEXP &get_data_from_simulation_SEXP,
                                           const List &simulation)
{
  GetDataFromSimulationPtr get_data_from_simulation = load_get_data_from_simulation(get_data_from_simulation_SEXP);
  return get_data_from_simulation(simulation);
}
