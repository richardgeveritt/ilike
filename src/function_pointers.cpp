#include <Rcpp.h>
using namespace Rcpp;

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
                                   const NumericVector &data)
{
  EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);
  return evaluate_log_likelihood(parameter, data);
}

// [[Rcpp::export]]
double estimate_log_likelihood_cpp(const SEXP &estimate_log_likelihood_SEXP,
                                   const NumericVector &parameter,
                                   const NumericVector &data,
                                   const List &auxiliary_variables)
{
  EstimateLogLikelihoodPtr estimate_log_likelihood = load_estimate_log_likelihood(estimate_log_likelihood_SEXP);
  return estimate_log_likelihood(parameter, data, auxiliary_variables);
}

// [[Rcpp::export]]
List simulate_model_cpp(const SEXP &simulate_model_SEXP,
                        const NumericVector &parameter,
                        const NumericVector &data)
{
  SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);
  return simulate_model(parameter, data);
}

// [[Rcpp::export]]
List simulate_auxiliary_variables_cpp(const SEXP &simulate_auxiliary_variables_SEXP,
                                      const NumericVector &parameter,
                                      const NumericVector &data)
{
  SimulateModelPtr simulate_auxiliary_variables = load_simulate_auxiliary_variables(simulate_auxiliary_variables_SEXP);
  return simulate_auxiliary_variables(parameter,data);
}
