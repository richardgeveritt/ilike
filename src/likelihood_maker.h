#include <Rcpp.h>
using namespace Rcpp;

#include "likelihood_estimator.h"

#ifndef LIKELIHOOD_MAKER
#define LIKELIHOOD_MAKER

LikelihoodEstimator* make_likelihood_estimator(const List &model,
                                               const List &algorithm);

#endif
