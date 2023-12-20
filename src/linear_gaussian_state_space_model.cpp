#include "linear_gaussian_state_space_model.h"
#include "utils.h"

/***data***/
Data data()
{
  arma::colvec Y(3);
  Y[0] = 0.0;
  Y[1] = 1.5;
  Y[2] = 1.8;
  
  Data data;
  data["y"] = Y;
  return data;
}

arma::colvec initial_mean()
{
  arma::colvec mean(2);
  mean[0] = 0.0;
  mean[1] = 0.0;
  return mean;
}

arma::mat initial_covariance()
{
  arma::mat covariance(2,2);
  covariance(0,0) = 10.0;
  covariance(1,1) = 10.0;
  return covariance;
}

/***evaluate_log_prior***/
double evaluate_log_prior(const Parameters &parameters)
{
  arma::colvec X = parameters["x"];
  return dmvnorm(X,initial_mean(),initial_covariance());
}

/***simulate_prior***/
Parameters simulate_prior(RandomNumberGenerator &rng)
{
  Parameters output;
  output["x"] = rmvnorm(rng,initial_mean(),initial_covariance());
  return output;
}

arma::mat transition_model_A(const Parameters &parameters)
{
  arma::mat A(2,2);
  A(0,0) = 1.0;
  A(0,1) = parameters["dt"][0];
  A(1,1) = 1.0;
  return A;
}

arma::mat transition_model_Q(const Parameters &parameters)
{
  arma::mat Q(2,2);
  Q(0,0) = sqrt(parameters["dt"][0]);
  Q(1,1) = sqrt(parameters["dt"][0]);
  return Q;
}

arma::mat measurement_model_H(const Parameters &parameters)
{
  arma::mat H(1,2);
  H(0,0) = 1.0;
  H(0,1) = 0.0;
  return H;

}

arma::mat measurement_model_R(const Parameters &parameters)
{
  arma::mat R(1,1);
  R(0,0) = 1.0;
  return R;
}
