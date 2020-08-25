#include <Rcpp.h>
using namespace Rcpp;

#include "exact_log_likelihood_estimator.h"
#include "function_pointers.h"

LogLikelihoodEstimator* make_log_likelihood_estimator(const List &model,
                                                      const List &algorithm)
{
  std::string likelihood_method = Rcpp::as<std::string>(algorithm["likelihood_method"]);
  if (likelihood_method=="analytic")
  {
    SEXP evaluate_log_likelihood_SEXP = model["evaluate_log_likelihood"];
    EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);
    return new ExactLogLikelihoodEstimator(evaluate_log_likelihood);
  }

  return NULL;
}
