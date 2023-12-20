#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ilike_header.h"

Data data();

arma::colvec initial_mean();
arma::mat initial_covariance();

double evaluate_log_prior(const Parameters &inputs);
Parameters simulate_prior(RandomNumberGenerator &rng);

arma::mat transition_model_A(const Parameters &parameters);
arma::mat transition_model_Q(const Parameters &parameters);
arma::mat measurement_model_H(const Parameters &parameters);
arma::mat measurement_model_R(const Parameters &parameters);
