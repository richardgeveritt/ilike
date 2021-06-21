//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "parameters.h"
#include "likelihood_estimator_output.h"

#ifndef LIKELIHOODESTIMATOR_H
#define LIKELIHOODESTIMATOR_H

class LikelihoodEstimator
{

protected:

  const Data* observed_data;

public:

  LikelihoodEstimator(const Data* observed_data_in);

  virtual ~LikelihoodEstimator();

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;

  virtual LikelihoodEstimatorOutput* simulate(const Parameters &parameters) const=0;

  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
};

#endif
