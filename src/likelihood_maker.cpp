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

    List observed_data = model["observed_data"];

    return new ExactLikelihoodEstimator(observed_data, evaluate_log_likelihood);
  }
  else if (likelihood_method=="abc")
  {
    SEXP simulate_model_SEXP = model["simulate_model"];
    SimulateModelPtr simulate_model = load_simulate_model(simulate_model_SEXP);

    //SEXP get_data_from_simulation_SEXP = algorithm["get_data_from_simulation"];
    //GetDataFromSimulationPtr get_data_from_simulation = load_get_data_from_simulation(get_data_from_simulation_SEXP);

    unsigned int number_of_likelihood_particles = algorithm["number_of_likelihood_particles"];

    SEXP evaluate_log_abc_kernel_SEXP = algorithm["evaluate_log_abc_kernel"];
    EvaluateLogABCKernelPtr evaluate_log_abc_kernel = load_evaluate_log_abc_kernel(evaluate_log_abc_kernel_SEXP);

    SEXP summary_statistics_SEXP = algorithm["summary_statistics"];
    SummaryStatisticsPtr summary_statistics = load_summary_statistics(summary_statistics_SEXP);

    double abc_tolerance = algorithm["abc_tolerance"];

    double abc_desired_cess = algorithm["abc_desired_cess"];

    arma::colvec summary_statistics_scaling = algorithm["summary_statistics_scaling"];

    bool adapt_abc_tolerance_to_cess = algorithm["adapt_abc_tolerance_to_cess"];

    bool adapt_summary_statistics_scaling = algorithm["adapt_summary_statistics_scaling"];

    List observed_data = model["observed_data"];

    return new ABCLikelihoodEstimator(observed_data,
                                      simulate_model,
                                      number_of_likelihood_particles,
                                      evaluate_log_abc_kernel,
                                      summary_statistics,
                                      abc_tolerance,
                                      abc_desired_cess,
                                      summary_statistics_scaling,
                                      adapt_abc_tolerance_to_cess,
                                      adapt_summary_statistics_scaling);
  }

  throw Rcpp::exception("likelihood_method not set to a valid option.");
}
