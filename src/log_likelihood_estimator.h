#include <Rcpp.h>
using namespace Rcpp;

#include <vector>

#ifndef LOGLIKELIHOODESTIMATOR_H
#define LOGLIKELIHOODESTIMATOR_H

class LogLikelihoodEstimator
{
public:
  LogLikelihoodEstimator();
  virtual ~LogLikelihoodEstimator();
  virtual double estimate_log_likelihood(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const=0;
  virtual List simulate_auxiliary_variables(const NumericVector &inputs,
                                            const NumericVector &data) const=0;
  virtual void setup_likelihood_estimator(const NumericMatrix &all_points,
                                          const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
