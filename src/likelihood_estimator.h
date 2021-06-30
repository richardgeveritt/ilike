#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"
#include "model_and_algorithm.h"

class LikelihoodEstimatorOutput;

class LikelihoodEstimator
{

public:

  LikelihoodEstimator();
  LikelihoodEstimator(const ModelAndAlgorithm &model_and_algorithm_in,
                      const Data* observed_data_in);
  virtual ~LikelihoodEstimator();

  LikelihoodEstimator(const LikelihoodEstimator &another);

  void operator=(const LikelihoodEstimator &another);
  virtual LikelihoodEstimator* duplicate() const=0;

  virtual LikelihoodEstimatorOutput* simulate(const Parameters &parameters)=0;

protected:

  // Not stored here.
  const Data* observed_data;

  // Stored here.
  ModelAndAlgorithm model_and_algorithm;

  void make_copy(const LikelihoodEstimator &another);

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;

  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
