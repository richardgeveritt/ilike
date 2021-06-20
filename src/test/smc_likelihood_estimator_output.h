#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "parameters.h"
#include "smc_output.h"

#ifndef SMCLIKELIHOODESTIMATOROUTPUT_H
#define SMCLIKELIHOODESTIMATOROUTPUT_H

class SMCLikelihoodEstimatorOutput
{

public:

  SMCLikelihoodEstimatorOutput();

  virtual ~SMCLikelihoodEstimatorOutput();

  void simulate(const Parameters &parameters);

  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  //
  // virtual List simulate_auxiliary_variables(const List &inputs) const=0;
  //
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;

protected:

  SMCOutput smc_output;

};

#endif

