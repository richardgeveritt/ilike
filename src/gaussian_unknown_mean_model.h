#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ilike_header.h"

Data data();

double evaluate_log_prior(const Parameters &inputs);

arma::mat evaluate_gradient_log_prior(const std::string &variable,
                                   const Parameters &parameters);

Parameters simulate_prior(RandomNumberGenerator &rng);

double evaluate_log_likelihood(const Parameters &inputs, const Data &observed_data);

arma::mat evaluate_gradient_log_likelihood(const std::string &variable,
                                           const Parameters &parameters,
                                           const Data &data);

Data simulate_model(RandomNumberGenerator &rng,
                    const Parameters &parameters);

Parameters simulate_independent_proposal(RandomNumberGenerator &rng);

double evaluate_log_proposal(const Parameters &inputs);

//double evaluate_log_abc_kernel(const Data &simulated_data, const Data &observed_data);

/*
Parameters transform_parameters(const Parameters &input);

Parameters inverse_transform_parameters(const Parameters &input);
*/

Data summary_statistics(const Data &data);
