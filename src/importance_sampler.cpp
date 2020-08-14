#include <Rcpp.h>
using namespace Rcpp;

#include "function_pointers.h"
#include "utils.h"


// [[Rcpp::export]]
List do_importance_sampler_cpp(const unsigned int &number_of_points,
                                        const List &model,
                                        List &algorithm,
                                        const unsigned int &max_vector_size)
{
  // Work out the number of batches to simulate.
  // We will store the points and the auxiliary variables in files to avoid memory problems.
  NumericVector inputs = model["inputs"];
  IntegerVector parameter_index = model["parameter_index"];
  NumericVector data = model["data"];
  unsigned int inputs_dimension = inputs.length();
  //unsigned int parameter_dimension = model["ParameterDimension"];
  unsigned int auxiliary_variables_dimension = algorithm["auxiliary_variables_dimension"];
  unsigned int total_dimension = inputs_dimension + auxiliary_variables_dimension;
  unsigned int number_of_batches_minus_one = floor(double(number_of_points*total_dimension-1)/max_vector_size);

  // Simulate batch of parameters and associated aux variables.
  // Save all in file if multiple batches.
  if (number_of_batches_minus_one==0)
  {
    // Everything can be done in one batch, so things are more straightforward.

    // Do the simulation.
    SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
    SimulateDistributionPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

    NumericMatrix proposed_inputs(number_of_points,inputs.size());
    NumericMatrix proposed_points(number_of_points,inputs.size());
    for (unsigned int i=0; i<number_of_points; ++i)
    {
      NumericVector proposed_point = simulate_proposal();
      NumericVector proposed_inputs_row = inputs;
      proposed_inputs_row[parameter_index] = proposed_point;
      proposed_points(i,_) = proposed_point;
      proposed_inputs(i,_) = proposed_inputs_row;
    }

    List likelihood_estimator = algorithm["likelihood_estimator"];

    bool evaluate_log_likelihood_is_set = likelihood_estimator["evaluate_log_likelihood_is_set"];

    List proposed_auxiliary_variables;

    if (evaluate_log_likelihood_is_set==FALSE)
    {
      SEXP simulate_auxiliary_variables_SEXP = likelihood_estimator["simulate_auxiliary_variables"];
      SimulateAuxiliaryVariablesPtr simulate_auxiliary_variables = load_simulate_auxiliary_variables(simulate_auxiliary_variables_SEXP);

      for (unsigned int i=0; i<number_of_points; ++i)
      {
        proposed_auxiliary_variables.push_back(simulate_auxiliary_variables(proposed_inputs(i,_),data));
      }

      // Now configure the likelihood estimator using all of the simulations, if needed.
      SEXP setup_likelihood_estimator_SEXP = likelihood_estimator["setup_likelihood_estimator"];
      SetupLikelihoodEstimatorPtr setup_likelihood_estimator = load_setup_likelihood_estimator(setup_likelihood_estimator_SEXP);

      XPtr<EvaluateLogLikelihoodPtr> evaluate_log_likelihood = setup_likelihood_estimator(proposed_points,proposed_auxiliary_variables);
      SEXP evaluate_log_likelihood_SEXP = SEXP(evaluate_log_likelihood);
      likelihood_estimator["evaluate_log_likelihood"] = evaluate_log_likelihood_SEXP;
      algorithm["likelihood_estimator"] = likelihood_estimator;
    }

    // Calculate weights.
    SEXP evaluate_log_likelihood_SEXP = likelihood_estimator["evaluate_log_likelihood"];
    EvaluateLogLikelihoodPtr evaluate_log_likelihood = load_evaluate_log_likelihood(evaluate_log_likelihood_SEXP);

    NumericVector log_weights(number_of_points);
    bool prior_is_proposal = algorithm["prior_is_proposal"];
    if (prior_is_proposal==TRUE)
    {
      for (unsigned int i=0; i<number_of_points; ++i)
      {
        log_weights[i] = evaluate_log_likelihood(proposed_inputs(i,_),data);
      }
    }
    else
    {
      SEXP evaluate_log_prior_SEXP = model["evaluate_log_prior"];
      EvaluateLogDistributionPtr evaluate_log_prior = load_evaluate_log_distribution(evaluate_log_prior_SEXP);

      SEXP evaluate_log_proposal_SEXP = algorithm["evaluate_log_proposal"];
      EvaluateLogDistributionPtr evaluate_log_proposal = load_evaluate_log_distribution(evaluate_log_proposal_SEXP);

      for (unsigned int i=0; i<number_of_points; ++i)
      {
        log_weights[i] = evaluate_log_likelihood(NumericVector(proposed_inputs(i,_)),data) + evaluate_log_prior(NumericVector(proposed_points(i,_))) - evaluate_log_proposal(NumericVector(proposed_points(i,_)));
      }
    }

    return List::create(Named("proposed_points") = proposed_points,
                        Named("proposed_auxiliary_variables") = proposed_auxiliary_variables,
                        Named("log_weights") = log_weights,
                        Named("log_normalising_constant") = log_sum_exp(log_weights));
  }
  else
  {
    throw Rcpp::exception("Cannot yet deal with multiple batches.");
  }

}
