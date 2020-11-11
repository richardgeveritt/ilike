#include "smc.h"

#include "function_pointers.h"
#include "utils.h"
#include "likelihood_estimator.h"
#include "likelihood_maker.h"

SMC::SMC(void)
{
}

SMC::SMC(const SMC &another)
{
  this->make_copy(another);
}

SMC::~SMC(void)
{
}

void SMC::operator=(const SMC &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void SMC::make_copy(const SMC &another)
{
  //Does nothing since no member variables to copy.
}

// [[Rcpp::export]]
List do_smc(const List &model,
            const List &algorithm)
{
  unsigned int number_of_points = algorithm["number_of_points"];

  List observed_data = model["observed_data"];

  // Do the initial importance sampling step.

  // Do the simulation.
  SEXP simulate_proposal_SEXP = algorithm["simulate_proposal"];
  SimulateDistributionPtr simulate_proposal = load_simulate_distribution(simulate_proposal_SEXP);

  std::vector<List> proposed_points;
  proposed_points.reserve(number_of_points);
  for (unsigned int i=0; i<number_of_points; ++i)
  {
    proposed_points.push_back(simulate_proposal());
  }

  LikelihoodEstimator* likelihood_estimator = make_likelihood_estimator(model, algorithm);

  std::vector<List> proposed_auxiliary_variables;
  proposed_auxiliary_variables.reserve(number_of_points);

  for (std::vector<List>::const_iterator i=proposed_points.begin(); i!=proposed_points.end(); ++i)
  {
    proposed_auxiliary_variables.push_back(likelihood_estimator->simulate_auxiliary_variables(*i));
  }

  likelihood_estimator->is_setup_likelihood_estimator(proposed_points,
                                                      proposed_auxiliary_variables);

  arma::colvec log_weights(number_of_points);
  bool prior_is_proposal = algorithm["prior_is_proposal"];
  if (prior_is_proposal==TRUE)
  {
    for (unsigned int i=0; i<number_of_points; ++i)
    {
      log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]);
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
      log_weights[i] = likelihood_estimator->estimate_log_likelihood(proposed_points[i], proposed_auxiliary_variables[i]) + evaluate_log_prior(proposed_points[i]) - evaluate_log_proposal(proposed_points[i]);
    }
  }

  if (likelihood_estimator != NULL)
    delete likelihood_estimator;

  return List::create(Named("proposed_points") = proposed_points,
                      Named("proposed_auxiliary_variables") = wrap(proposed_auxiliary_variables),
                      Named("log_weights") = log_weights,
                      Named("log_normalising_constant") = log_sum_exp(log_weights));
}
