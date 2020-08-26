#include <Rcpp.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "abc_likelihood.h"

#ifndef ABCLIKELIHOODESTIMATOR_H
#define ABCLIKELIHOODESTIMATOR_H

class ABCLikelihoodEstimator : public LikelihoodEstimator
{

private:

  SimulateModelPtr simulator;

  unsigned int number_of_simulations;

  ABCLikelihood abc_likelihood;

public:

  ABCLikelihoodEstimator(const NumericVector &data_in,
                         const SimulateModelPtr &simulator_in,
                         const double &number_of_simulations_in,
                         const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                         const SummaryStatisticPtr &summary_statistic_in,
                         const NumericVector &abc_tolerances_in,
                         const NumericVector &summary_statistic_scaling_in);

  virtual ~ABCLikelihoodEstimator();

  double estimate_log_likelihood(const NumericVector &inputs,
                                 const List &auxiliary_variables) const;

  List simulate_auxiliary_variables(const NumericVector &inputs) const;

  void setup_likelihood_estimator(const NumericMatrix &all_points,
                                  const std::vector<List> &all_auxiliary_variables);
};

#endif
