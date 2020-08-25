#include <Rcpp.h>
using namespace Rcpp;

#include "log_likelihood_estimator.h"

#ifndef LIKELIHOOD_MAKER
#define LIKELIHOOD_MAKER

LogLikelihoodEstimator* make_log_likelihood_estimator(const List &model,
                                                      const List &algorithm);

#endif
