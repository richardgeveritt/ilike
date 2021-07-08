#ifndef LIKELIHOODESTIMATOROUTPUT_H
#define LIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

class LikelihoodEstimatorOutput
{

public:

  LikelihoodEstimatorOutput();

  virtual ~LikelihoodEstimatorOutput();

  //virtual void simulate(const Parameters &parameters);

  LikelihoodEstimatorOutput(const LikelihoodEstimatorOutput &another);

  void operator=(const LikelihoodEstimatorOutput &another);
  virtual LikelihoodEstimatorOutput* duplicate() const=0;

  virtual void continue_simulate(const Parameters &parameters)=0;
  virtual void estimate(const Parameters &parameters)=0;

  double log_likelihood;

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  //
  // virtual List simulate_auxiliary_variables(const List &inputs) const=0;
  //
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;

protected:

  void make_copy(const LikelihoodEstimatorOutput &another);

};

#endif

