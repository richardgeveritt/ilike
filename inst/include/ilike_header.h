#ifndef TYPEDEFS
#define TYPEDEFS

// [[Rcpp::interfaces(r, cpp)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"
#include "distributions.h"
//#include "data.h"

// typedefs for function pointers

namespace ilike
{

typedef double (*EvaluateLogDistributionPtr)(const Parameters &parameters);

typedef Parameters (*SimulateDistributionPtr)(RandomNumberGenerator &rng);

typedef double (*EvaluateLogLikelihoodPtr)(const Parameters &parameters,
                                           const Data &data);

typedef arma::mat (*EvaluateGradientLogDistributionPtr)(const std::string &variable,
                                                        const Parameters &parameters);

typedef arma::mat (*EvaluateGradientLogLikelihoodPtr)(const std::string &variable,
                                                      const Parameters &parameters,
                                                      const Data &data);

typedef double (*EvaluateLogNoParamsMCMCProposalPtr)(const Parameters &proposed_parameters,
                                                     const Parameters &parameters);

typedef Parameters (*SimulateNoParamsMCMCProposalPtr)(RandomNumberGenerator &rng,
                                                      const Parameters &parameters);

typedef double (*EvaluateLogMCMCProposalPtr)(const Parameters &proposed_parameters,
                                             const Parameters &parameters,
                                             const Parameters &proposal_parameters);

typedef Parameters (*SimulateMCMCProposalPtr)(RandomNumberGenerator &rng,
                                              const Parameters &parameters,
                                              const Parameters &proposal_parameters);

typedef double (*EvaluateLogGuidedNoParamsMCMCProposalPtr)(const Parameters &proposed_parameters,
                                                           const Parameters &parameters,
                                                           const Data &data);

typedef Parameters (*SimulateGuidedNoParamsMCMCProposalPtr)(RandomNumberGenerator &rng,
                                                            const Parameters &parameters,
                                                            const Data &data);

typedef double (*EvaluateLogGuidedMCMCProposalPtr)(const Parameters &proposed_parameters,
                                                   const Parameters &parameters,
                                                   const Parameters &proposal_parameters,
                                                   const Data &data);

typedef Parameters (*SimulateGuidedMCMCProposalPtr)(RandomNumberGenerator &rng,
                                                    const Parameters &parameters,
                                                    const Parameters &proposal_parameters,
                                                    const Data &data);

typedef double (*EvaluateLogIndependentProposalPtr)(const Parameters &proposed_parameters,
                                                    const Parameters &proposal_parameters);

typedef Parameters (*SimulateIndependentProposalPtr)(RandomNumberGenerator &rng,
                                                     const Parameters &proposal_parameters);

typedef double (*EvaluateLogGuidedDistributionPtr)(const Parameters &proposed_parameters,
                                                   const Data &data);

typedef Parameters (*SimulateGuidedDistributionPtr)(RandomNumberGenerator &rng,
                                                    const Data &data);

typedef double (*EvaluateLogIndependentProposalPtr)(const Parameters &proposed_parameters,
                                                    const Parameters &proposal_parameters);

typedef Parameters (*SimulateIndependentProposalPtr)(RandomNumberGenerator &rng,
                                                     const Parameters &proposal_parameters);

typedef double (*EvaluateLogGuidedIndependentProposalPtr)(const Parameters &proposed_parameters,
                                                          const Parameters &proposal_parameters,
                                                          const Data &data);

typedef Parameters (*SimulateGuidedIndependentProposalPtr)(RandomNumberGenerator &rng,
                                                           const Parameters &proposal_parameters,
                                                           const Data &data);

typedef Data (*SimulateModelPtr)(RandomNumberGenerator &rng,
                                 const Parameters &parameters);

typedef double (*GainPtr)(size_t n);

typedef Parameters (*TransformPtr)(const Parameters &inputs);

typedef Data (*TransformToDataPtr)(const Parameters &inputs);

typedef Data (*SummaryStatisticsPtr)(const Data &inputs);

typedef arma::mat (*JacobianPtr)(const Parameters &inputs);

typedef arma::mat (*GetProcessMatrixFromTimePtr)(double time_step);
//typedef arma::mat (*GetProcessMatrixFromParametersPtr)(const Parameters &conditioned_on_parameters);
typedef arma::mat (*GetProcessMatrixFromTimeParametersPtr)(double time_step,
                                                           const Parameters &conditioned_on_parameters);

typedef arma::colvec (*SimulateTransitionKernelPtr)(const arma::colvec &current_state);
typedef arma::colvec (*SimulateTransitionKernelFromTimePtr)(const arma::colvec &current_state,
                                                            double time_step);
typedef arma::colvec (*SimulateTransitionKernelFromParametersPtr)(const arma::colvec &current_state,
                                                                  const Parameters &conditioned_on_parameters);
typedef arma::colvec (*SimulateTransitionKernelFromTimeParametersPtr)(const arma::colvec &current_state,
                                                                      double time_step,
                                                                      const Parameters &conditioned_on_parameters);

typedef arma::colvec (*GetNoParamsVectorPtr)();

typedef arma::mat (*GetNoParamsMatrixPtr)();

typedef arma::colvec (*GetVectorPtr)(const Parameters &conditioned_on_parameters);

typedef arma::mat (*GetMatrixPtr)(const Parameters &conditioned_on_parameters);

typedef arma::colvec (*MatrixSimulateMeasurementKernelPtr)(const arma::colvec &current_state);

typedef arma::colvec (*GetMatrixSimulateMeasurementKernelPtr)(const arma::colvec &current_state,
                                                              const Parameters &conditioned_on_parameters);

typedef Parameters (*SimulateMeasurementKernelPtr)(const Parameters &current_state);

typedef Parameters (*GetSimulateMeasurementKernelPtr)(const Parameters &current_state,
                                                      const Parameters &conditioned_on_parameters);

typedef double (*PowerFunctionPtr)(const Parameters &current_state,
                                   const std::string &power_variable);



typedef double (*EstimateLogLikelihoodPtr)(const List &inputs, const List &observed_data, const List &auxiliary_variables);

//typedef List (*SimulateModelPtr)(const List &inputs, const List &observed_data);

typedef List (*SimulateAuxiliaryVariablesPtr)(const List &inputs, const List &observed_data);

typedef XPtr<EstimateLogLikelihoodPtr> (*SetupLikelihoodEstimatorPtr)(const List &inputs, const List &auxiliary_variables);

typedef double (*EvaluateLogABCKernelPtr)(const arma::colvec &simulated_stats,
                                          const arma::colvec &observed_stats,
                                          const double &abc_tolerance);

//typedef arma::colvec (*SummaryStatisticsPtr)(const List &observed_data);

typedef List (*GetDataFromSimulationPtr)(const List &simulation);

typedef Data (*DataPtr)();

typedef NumericVector (*NumericVectorPtr)();

typedef std::string (*StringPtr)();


/*
 EvaluateLogDistributionPtr load_evaluate_log_distribution(const SEXP &evaluate_log_distribution_SEXP);
 
 SimulateDistributionPtr load_simulate_distribution(const SEXP &simulate_distribution_SEXP);
 
 EvaluateLogGuidedDistributionPtr load_evaluate_log_guided_distribution(const SEXP &evaluate_log_guided_distribution_SEXP);
 
 SimulateGuidedDistributionPtr load_simulate_guided_distribution(const SEXP &simulate_guided_distribution_SEXP);
 
 EvaluateLogIndependentProposalPtr load_evaluate_log_independent_proposal(const SEXP &evaluate_log_independent_proposal_SEXP);
 
 SimulateIndependentProposalPtr load_simulate_independent_proposal(const SEXP &simulate_independent_proposal_SEXP);
 
 EvaluateLogGuidedIndependentProposalPtr load_evaluate_log_guided_independent_proposal(const SEXP &evaluate_log_guided_independent_proposal_SEXP);
 
 SimulateGuidedIndependentProposalPtr load_simulate_guided_independent_proposal(const SEXP &simulate_guided_independent_proposal_SEXP);
 
 EvaluateLogNoParamsMCMCProposalPtr load_evaluate_log_no_params_mcmc_proposal(const SEXP &evaluate_log_no_params_mcmc_proposal_SEXP);
 
 SimulateNoParamsMCMCProposalPtr load_simulate_no_params_mcmc_proposal(const SEXP &simulate_no_params_mcmc_proposal_SEXP);
 
 EvaluateLogGuidedNoParamsMCMCProposalPtr load_evaluate_log_guided_no_params_mcmc_proposal(const SEXP &evaluate_log_guided_no_params_mcmc_proposal_SEXP);
 
 SimulateGuidedNoParamsMCMCProposalPtr load_simulate_guided_no_params_mcmc_proposal(const SEXP &simulate_guided_no_params_mcmc_proposal_SEXP);
 
 EvaluateLogMCMCProposalPtr load_evaluate_log_mcmc_proposal(const SEXP &evaluate_log_mcmc_proposal_SEXP);
 
 SimulateMCMCProposalPtr load_simulate_mcmc_proposal(const SEXP &simulate_mcmc_proposal_SEXP);
 
 EvaluateLogGuidedMCMCProposalPtr load_evaluate_log_guided_mcmc_proposal(const SEXP &evaluate_log_guided_mcmc_proposal_SEXP);
 
 SimulateGuidedMCMCProposalPtr load_simulate_guided_mcmc_proposal(const SEXP &simulate_mcmc_guided_proposal_SEXP);
 
 EvaluateLogLikelihoodPtr load_evaluate_log_likelihood(const SEXP &evaluate_log_likelihood_SEXP);
 
 SimulateModelPtr load_simulate_data_model(const SEXP &simulate_data_model_SEXP);
 
 SimulateAuxiliaryVariablesPtr load_simulate_auxiliary_variables(const SEXP &simulate_auxiliary_variables_SEXP);
 
 SetupLikelihoodEstimatorPtr load_setup_likelihood_estimator(const SEXP &setup_likelihood_estimator_SEXP);
 
 EvaluateLogABCKernelPtr load_evaluate_log_abc_kernel(const SEXP &evaluate_log_abc_kernel_SEXP);
 
 SummaryStatisticsPtr load_summary_statistics(const SEXP &summary_statistic_SEXP);
 
 GetDataFromSimulationPtr load_get_data_from_simulation(const SEXP &get_data_from_simulation_SEXP);
 
 //Data load_observed_data(const SEXP &observed_data_SEXP);
 
 NumericVector load_mcmc_weights(const SEXP &mcmc_weights_SEXP);
 
 Data load_data(const SEXP &data_SEXP);
 
 DataPtr load_data_function(const SEXP &data_SEXP);
 */

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
// List simulate_data_model_cpp(const SEXP &simulate_data_model_SEXP,
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

inline EvaluateLogGuidedDistributionPtr load_evaluate_log_guided_distribution(const SEXP &evaluate_log_guided_distribution_SEXP)
{
  XPtr<EvaluateLogGuidedDistributionPtr> evaluate_log_guided_distribution_XPtr(evaluate_log_guided_distribution_SEXP);
  return *evaluate_log_guided_distribution_XPtr;
}

inline SimulateGuidedDistributionPtr load_simulate_guided_distribution(const SEXP &simulate_guided_distribution_SEXP)
{
  XPtr<SimulateGuidedDistributionPtr> simulate_guided_distribution_XPtr(simulate_guided_distribution_SEXP);
  return *simulate_guided_distribution_XPtr;
}

inline EvaluateLogIndependentProposalPtr load_evaluate_log_independent_proposal(const SEXP &evaluate_log_independent_proposal_SEXP)
{
  XPtr<EvaluateLogIndependentProposalPtr> evaluate_log_independent_proposal_XPtr(evaluate_log_independent_proposal_SEXP);
  return *evaluate_log_independent_proposal_XPtr;
}

inline SimulateGuidedIndependentProposalPtr load_guided_simulate_independent_proposal(const SEXP &simulate_guided_independent_proposal_SEXP)
{
  XPtr<SimulateGuidedIndependentProposalPtr> simulate_guided_distribution_XPtr(simulate_guided_independent_proposal_SEXP);
  return *simulate_guided_distribution_XPtr;
}

inline EvaluateLogGuidedIndependentProposalPtr load_evaluate_log_guided_independent_proposal(const SEXP &evaluate_log_guided_independent_proposal_SEXP)
{
  XPtr<EvaluateLogGuidedIndependentProposalPtr> evaluate_log_guided_independent_proposal_XPtr(evaluate_log_guided_independent_proposal_SEXP);
  return *evaluate_log_guided_independent_proposal_XPtr;
}

inline SimulateGuidedIndependentProposalPtr load_simulate_guided_independent_proposal(const SEXP &simulate_guided_independent_proposal_SEXP)
{
  XPtr<SimulateGuidedIndependentProposalPtr> simulate_guided_independent_proposal_XPtr(simulate_guided_independent_proposal_SEXP);
  return *simulate_guided_independent_proposal_XPtr;
}

inline SimulateIndependentProposalPtr load_simulate_independent_proposal(const SEXP &simulate_independent_proposal_SEXP)
{
  XPtr<SimulateIndependentProposalPtr> simulate_independent_proposal_XPtr(simulate_independent_proposal_SEXP);
  return *simulate_independent_proposal_XPtr;
}

inline EvaluateLogNoParamsMCMCProposalPtr load_evaluate_log_no_params_mcmc_proposal(const SEXP &evaluate_log_no_params_mcmc_proposal_SEXP)
{
  XPtr<EvaluateLogNoParamsMCMCProposalPtr> evaluate_log_no_params_mcmc_proposal_XPtr(evaluate_log_no_params_mcmc_proposal_SEXP);
  return *evaluate_log_no_params_mcmc_proposal_XPtr;
}

inline SimulateNoParamsMCMCProposalPtr load_simulate_no_params_mcmc_proposal(const SEXP &simulate_no_params_mcmc_proposal_SEXP)
{
  XPtr<SimulateNoParamsMCMCProposalPtr> simulate_no_params_mcmc_proposal_XPtr(simulate_no_params_mcmc_proposal_SEXP);
  return *simulate_no_params_mcmc_proposal_XPtr;
}

inline EvaluateLogGuidedNoParamsMCMCProposalPtr load_evaluate_log_guided_no_params_mcmc_proposal(const SEXP &evaluate_log_mcmc_guided_no_params_proposal_SEXP)
{
  XPtr<EvaluateLogGuidedNoParamsMCMCProposalPtr> evaluate_log_mcmc_guided_no_params_proposal_XPtr(evaluate_log_mcmc_guided_no_params_proposal_SEXP);
  return *evaluate_log_mcmc_guided_no_params_proposal_XPtr;
}

inline SimulateGuidedNoParamsMCMCProposalPtr load_simulate_guided_no_params_mcmc_proposal(const SEXP &simulate_mcmc_guided_no_params_proposal_SEXP)
{
  XPtr<SimulateGuidedNoParamsMCMCProposalPtr> simulate_mcmc_guided_no_params_proposal_XPtr(simulate_mcmc_guided_no_params_proposal_SEXP);
  return *simulate_mcmc_guided_no_params_proposal_XPtr;
}

inline EvaluateLogMCMCProposalPtr load_evaluate_log_mcmc_proposal(const SEXP &evaluate_log_mcmc_proposal_SEXP)
{
  XPtr<EvaluateLogMCMCProposalPtr> evaluate_log_mcmc_proposal_XPtr(evaluate_log_mcmc_proposal_SEXP);
  return *evaluate_log_mcmc_proposal_XPtr;
}

inline SimulateMCMCProposalPtr load_simulate_mcmc_proposal(const SEXP &simulate_mcmc_proposal_SEXP)
{
  XPtr<SimulateMCMCProposalPtr> simulate_mcmc_proposal_XPtr(simulate_mcmc_proposal_SEXP);
  return *simulate_mcmc_proposal_XPtr;
}

inline EvaluateLogGuidedMCMCProposalPtr load_evaluate_log_guided_mcmc_proposal(const SEXP &evaluate_log_guided_mcmc_proposal_SEXP)
{
  XPtr<EvaluateLogGuidedMCMCProposalPtr> evaluate_log_guided_mcmc_proposal_XPtr(evaluate_log_guided_mcmc_proposal_SEXP);
  return *evaluate_log_guided_mcmc_proposal_XPtr;
}

inline SimulateGuidedMCMCProposalPtr load_simulate_guided_mcmc_proposal(const SEXP &simulate_guided_mcmc_proposal_SEXP)
{
  XPtr<SimulateGuidedMCMCProposalPtr> simulate_guided_mcmc_proposal_XPtr(simulate_guided_mcmc_proposal_SEXP);
  return *simulate_guided_mcmc_proposal_XPtr;
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

inline SimulateModelPtr load_simulate_data_model(const SEXP &simulate_data_model_SEXP)
{
  XPtr<SimulateModelPtr> simulate_data_model_XPtr(simulate_data_model_SEXP);
  return *simulate_data_model_XPtr;
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

inline TransformPtr load_transform(const SEXP &transform_SEXP)
{
  XPtr<TransformPtr> transform_XPtr(transform_SEXP);
  return *transform_XPtr;
}

inline GetDataFromSimulationPtr load_get_data_from_simulation(const SEXP &get_data_from_simulation_SEXP)
{
  XPtr<GetDataFromSimulationPtr> get_data_from_simulation_XPtr(get_data_from_simulation_SEXP);
  return *get_data_from_simulation_XPtr;
}

inline arma::colvec load_vector(const SEXP &get_vector_SEXP)
{
  XPtr<GetNoParamsVectorPtr> get_vector_XPtr(get_vector_SEXP);
  return (*get_vector_XPtr)();
}

inline arma::mat load_matrix(const SEXP &get_matrix_SEXP)
{
  XPtr<GetNoParamsMatrixPtr> get_matrix_XPtr(get_matrix_SEXP);
  return (*get_matrix_XPtr)();
}

inline GetVectorPtr load_vector_function(const SEXP &get_vector_SEXP)
{
  XPtr<GetVectorPtr> get_vector_XPtr(get_vector_SEXP);
  return *get_vector_XPtr;
}

inline GetMatrixPtr load_matrix_function(const SEXP &get_matrix_SEXP)
{
  XPtr<GetMatrixPtr> get_matrix_XPtr(get_matrix_SEXP);
  return *get_matrix_XPtr;
}

// inline Data load_observed_data(const SEXP &observed_data_SEXP)
// {
//   XPtr<Data> observed_data_XPtr(observed_data_SEXP);
//   return *observed_data_XPtr;
// }

inline NumericVector load_numeric_vector(const SEXP &numeric_vector_SEXP)
{
  XPtr<NumericVectorPtr> numeric_vector_XPtr(numeric_vector_SEXP);
  //return Parameters();
  return (*numeric_vector_XPtr)();
}

inline Data load_data(const SEXP &data_SEXP)
{
  XPtr<DataPtr> data_XPtr(data_SEXP);
  //return Parameters();
  return (*data_XPtr)();
}

inline DataPtr load_data_function(const SEXP &data_SEXP)
{
  XPtr<DataPtr> data_XPtr(data_SEXP);
  //return Parameters();
  return *data_XPtr;
}

inline std::string load_string(const SEXP &string_SEXP)
{
  XPtr<StringPtr> string_XPtr(string_SEXP);
  //return Parameters();
  return (*string_XPtr)();
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
// List simulate_data_model_cpp(const SEXP &simulate_data_model_SEXP,
//                         const List &parameter,
//                         const List &observed_data)
// {
//   SimulateModelPtr simulate_data_model = load_simulate_data_model(simulate_data_model_SEXP);
//   return simulate_data_model(parameter, observed_data);
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

}
#endif
