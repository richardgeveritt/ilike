#ifndef EXACTLIKELIHOODESTIMATOR_H
#define EXACTLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"

class ExactLikelihoodEstimatorOutput;

class ExactLikelihoodEstimator : public LikelihoodEstimator
{

public:

  ExactLikelihoodEstimator();

  ExactLikelihoodEstimator(RandomNumberGenerator* rng_in,
                           size_t* seed_in,
                           const Data* sdata_in,
                           EvaluateLogLikelihoodPtr func_in);

  virtual ~ExactLikelihoodEstimator();

  ExactLikelihoodEstimator(const ExactLikelihoodEstimator &another);

  void operator=(const ExactLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* initial_simulate(const Parameters &parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

private:

  friend ExactLikelihoodEstimatorOutput;

  EvaluateLogLikelihoodPtr func;

  // Stored here.
  //ExactLikelihoodEstimatorOutput* output;

  void make_copy(const ExactLikelihoodEstimator &another);

};

#endif
