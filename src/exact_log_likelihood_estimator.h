#include <Rcpp.h>
using namespace Rcpp;

#include <vector>

#include "log_likelihood_estimator.h"
#include "function_pointers.h"

#ifndef EXACTLOGLIKELIHOODESTIMATOR_H
#define EXACTLOGLIKELIHOODESTIMATOR_H

class ExactLogLikelihoodEstimator : public LogLikelihoodEstimator
{
private:
  EvaluateLogLikelihoodPtr func;
public:
  ExactLogLikelihoodEstimator(const EvaluateLogLikelihoodPtr &func_in);
  virtual ~ExactLogLikelihoodEstimator();
  double estimate_log_likelihood(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const;
  List simulate_auxiliary_variables(const NumericVector &inputs,
                                    const NumericVector &data) const;
  void setup_likelihood_estimator(const NumericMatrix &all_points,
                                  const std::vector<List> &all_auxiliary_variables);
};

#endif
