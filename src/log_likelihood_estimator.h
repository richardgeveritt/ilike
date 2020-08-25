#include <Rcpp.h>
using namespace Rcpp;

#ifndef LOGLIKELIHOODESTIMATOR_H
#define LOGLIKELIHOODESTIMATOR_H

class LogLikelihoodEstimator
{
public:
  LogLikelihoodEstimator();
  virtual ~LogLikelihoodEstimator();
  virtual double operator()(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const=0;
  virtual List simulate_auxiliary_variables() const=0;
};

#endif
