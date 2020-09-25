#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"

#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

class ExactLikelihoodEstimator : public LikelihoodEstimator
{

private:

  EvaluateLogLikelihoodPtr func;

public:

  ExactLikelihoodEstimator(const List &observed_data_in,
                           const EvaluateLogLikelihoodPtr &func_in);

  virtual ~ExactLikelihoodEstimator();

  double estimate_log_likelihood(const List &inputs,
                                 const List &auxiliary_variables) const;

  List simulate_auxiliary_variables(const List &inputs) const;

  void is_setup_likelihood_estimator(const std::vector<List> &all_points,
                                     const std::vector<List> &all_auxiliary_variables);
};

#endif
