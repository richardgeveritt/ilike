#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
#include "model_and_algorithm.h"
#include "distributions.h"

class LikelihoodEstimatorOutput;
class SMCWorker;

class LikelihoodEstimator
{

public:

  LikelihoodEstimator();
  LikelihoodEstimator(RandomNumberGenerator* rng_in,
                      size_t* seed_in,
                      const Data* data_in);
  virtual ~LikelihoodEstimator();

  LikelihoodEstimator(const LikelihoodEstimator &another);

  void operator=(const LikelihoodEstimator &another);
  virtual LikelihoodEstimator* duplicate() const=0;

  virtual LikelihoodEstimatorOutput* initial_simulate(const Parameters &parameters)=0;

protected:

  friend SMCWorker;

  // Not stored here. Stored in "main'.
  const Data* data;

  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;

  // Not stored here. Stored in "main'.
  size_t* seed;

  ModelAndAlgorithm model_and_algorithm;

  void make_copy(const LikelihoodEstimator &another);

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;

  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
