//#include <Rcpp.h>
//using namespace Rcpp;

#include <RcppArmadillo.h>
using namespace Rcpp;


// typedefs for function pointers

#ifndef TYPEDEFS
#define TYPEDEFS


typedef double (*EvaluateLogDistributionPtr)(const NumericVector &inputs);

typedef NumericVector (*SimulateDistributionPtr)(void);

typedef double (*EvaluateLogLikelihoodPtr)(const NumericVector &inputs, const NumericMatrix &data);

typedef double (*EstimateLogLikelihoodPtr)(const NumericVector &inputs, const NumericMatrix &data, const List &auxiliary_variables);

typedef List (*SimulateModelPtr)(const NumericVector &inputs, const NumericMatrix &data);

typedef List (*SimulateAuxiliaryVariablesPtr)(const NumericVector &inputs, const NumericMatrix &data);

typedef XPtr<EstimateLogLikelihoodPtr> (*SetupLikelihoodEstimatorPtr)(const NumericMatrix &inputs, const List &auxiliary_variables);

typedef double (*EvaluateLogABCKernelPtr)(const arma::colvec &simulated_stats,
                const arma::colvec &observed_stats,
                const double &abc_tolerance);

typedef arma::colvec (*SummaryStatisticsPtr)(const NumericMatrix &data);

typedef NumericMatrix (*GetDataFromSimulationPtr)(const List &simulation);


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


NumericVector simulate_distribution_cpp(const SEXP &simulate_distribution_SEXP);

double evaluate_log_distribution_cpp(const SEXP &evaluate_log_distribution_SEXP,
                                     const NumericVector &parameter);

double evaluate_log_likelihood_cpp(const SEXP &evaluate_log_likelihood_SEXP,
                                   const NumericVector &parameter,
                                   const NumericMatrix &data);

double estimate_log_likelihood_cpp(const SEXP &estimate_log_likelihood_SEXP,
                                   const NumericVector &parameter,
                                   const NumericMatrix &data,
                                   const List &auxiliary_variables);

List simulate_model_cpp(const SEXP &simulate_model_SEXP,
                        const NumericVector &parameter,
                        const NumericMatrix &data);

List simulate_auxiliary_variables_cpp(const SEXP &simulate_auxiliary_variables_SEXP,
                                      const NumericVector &parameter,
                                      const NumericMatrix &data);

double evaluate_log_abc_kernel_cpp(const SEXP &evaluate_log_abc_kernel_SEXP,
                                   const arma::colvec &simulated_data,
                                   const arma::colvec &data,
                                   const double &abc_tolerance);

arma::colvec summary_statistics_cpp(const SEXP &summary_statistics_SEXP,
                                    const NumericMatrix &data);

NumericMatrix get_data_from_simulation_cpp(const SEXP &get_data_from_simulation_SEXP,
                                           const List &simulation);

#endif
