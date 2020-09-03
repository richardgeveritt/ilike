#include <Rcpp.h>
using namespace Rcpp;

#include "exact_likelihood_estimator.h"
#include "abc_likelihood_estimator.h"
#include "function_pointers.h"

LikelihoodEstimator* make_likelihood_estimator(const List &model,
                                               const List &algorithm)
{
  std::string likelihood_method = Rcpp::as<std::string>(algorithm["likelihood_method"]);
  if (likelihood_method=="analytic")
  {
    SEXP evaluate_log_likelihood_SEXP = model["evaluate_log_likelihood"];
    EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);

    NumericVector data = model["data"];

    return new ExactLikelihoodEstimator(data, evaluate_log_likelihood);
  }
  else if (likelihood_method=="abc")
  {
    SEXP simulate_model_SEXP = model["simulate_model"];
    SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);

    unsigned int number_of_simulations = algorithm["number_of_simulations"];

    SEXP evaluate_log_abc_kernel_SEXP = algorithm["evaluate_log_abc_kernel"];
    EvaluateLogABCKernelPtr evaluate_log_abc_kernel = load_evaluate_log_abc_kernel(evaluate_log_abc_kernel_SEXP);

    SEXP summary_statistics_SEXP = algorithm["summary_statistics"];
    SummaryStatisticsPtr summary_statistics = load_summary_statistics(summary_statistics_SEXP);

    double abc_tolerance = algorithm["abc_tolerance"];

    NumericVector summary_statistics_scaling = algorithm["summary_statistics_scaling"];

    NumericVector data = model["data"];

    return new ABCLikelihoodEstimator(data,
                                      simulate_model,
                                      number_of_simulations,
                                      evaluate_log_abc_kernel,
                                      summary_statistics,
                                      abc_tolerance,
                                      summary_statistics_scaling);
  }

  return NULL;
}
