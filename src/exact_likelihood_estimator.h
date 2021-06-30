#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"

class LikelihoodEstimatorOutput;

class ExactLikelihoodEstimator : public LikelihoodEstimator
{

public:

  ExactLikelihoodEstimator(const ModelAndAlgorithm &model_and_algorithm_in,
                           const Data* observed_data_in,
                           EvaluateLogLikelihoodPtr func_in);

  virtual ~ExactLikelihoodEstimator();

  ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another);

  void operator=(const ExactLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* simulate(const Parameters &parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

private:

  EvaluateLogLikelihoodPtr func;

  void make_copy(const ExactLikelihoodEstimator &another);

};

#endif
