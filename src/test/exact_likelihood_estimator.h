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

  ExactLikelihoodEstimator(const Data* observed_data_in,
                           const EvaluateLogLikelihoodPtr &func_in);

  virtual ~ExactLikelihoodEstimator();

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* simulate(const Parameters &parameters) const;

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
};

#endif
