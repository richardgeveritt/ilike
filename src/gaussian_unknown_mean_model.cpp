#include "gaussian_unknown_mean_model.h"
#include "utils.h"

/***data***/
Data data()
{
  RandomNumberGenerator rng;
  rng.seed(100);
  size_t n = 100;
  arma::colvec sampled = rnorm(rng,n);
  Data data;
  data["y"] = sampled;
  return data;
}

/***evaluate_log_prior***/
double evaluate_log_prior(const Parameters &parameters)
{
  return dnorm(parameters["mu"][0], 0.0, 10.0);
}

/***evaluate_gradient_log_prior***/
arma::mat evaluate_gradient_log_prior(const std::string &variable,
                                      const Parameters &parameters)
{
  if (variable=="mu")
  {
    return -parameters[variable];
  }
  else
  {
    stop("Can only take gradient with respect to mu.");
  }
}

/***simulate_prior***/
Parameters simulate_prior(RandomNumberGenerator &rng)
{
  Parameters output;
  output["mu"] = rnorm(rng, 0.0, 10.0);
  return output;
}

/***evaluate_log_likelihood***/
double evaluate_log_likelihood(const Parameters &parameters, const Data &data)
{
  double mean = parameters["mu"][0];
  return sum(dnorm(data["y"], mean, 1.0));
}

/***evaluate_gradient_log_likelihood***/
arma::mat evaluate_gradient_log_likelihood(const std::string &variable,
                                           const Parameters &parameters,
                                           const Data &data)
{
  if (variable=="mu")
  {
    return -arma::sum(data["y"]-parameters[variable][0]);
  }
  else
  {
    stop("Can only take gradient with respect to mu.");
  }
}

/***simulate_model***/
Data simulate_model(RandomNumberGenerator &rng,
                    const Parameters &parameters)
{
  Data output;
  size_t n = 100;
  double mean = parameters["mean"][0];
  output["y"] = rnorm(rng, n, mean, 1.0);
  return output;
}

/***simulate_importance_proposal***/
Parameters simulate_importance_proposal(RandomNumberGenerator &rng,
                                        const Parameters &proposal_parameters)
{
  Parameters output;
  output["mean"] = rnorm(rng, 0.0, 10.0);
  return output;
}

/***evaluate_log_importance_proposal***/
double evaluate_log_importance_proposal(const Parameters &parameters)
{
  return dnorm(parameters["mean"][0], 0.0, 10.0);
}

/*
double evaluate_log_abc_kernel(const Data &simulated_data, const Data &observed_data)
{
  return uniform_abc_kernel(simulated_data["y"],observed_data["y"],simulated_data["epsilon"][0]);
}
*/

/***summary_statistics***/
Data summary_statistics(const Data &data)
{
  return Data("mean",arma::mean(data["y"]));
}
