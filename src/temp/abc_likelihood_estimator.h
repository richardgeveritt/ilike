#include <RcppArmadillo.h>
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

  unsigned int number_of_likelihood_particles;

  ABCLikelihood abc_likelihood;

  double abc_desired_cess;

  bool adapt_abc_tolerance_to_cess;

  bool adapt_summary_statistics_scaling;

  double cess_score(const double &current_epsilon,
                    const std::vector<List> &all_points,
                    const std::vector<List> &all_auxiliary_variables,
                    const arma::colvec &current_log_weights);

  double epsilon_bisection(const double &current_epsilon,
                           const std::vector<List> &all_points,
                           const std::vector<List> &all_auxiliary_variables,
                           const arma::colvec &current_log_weights);

  double epsilon_doubling(const std::vector<List> &all_points,
                          const std::vector<List> &all_auxiliary_variables,
                          const arma::colvec &current_log_weights);

public:

  ABCLikelihoodEstimator(const List &observed_data_in,
                         const SimulateModelPtr &simulator_in,
                         const double &number_of_likelihood_particles_in,
                         const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                         const SummaryStatisticsPtr &summary_statistics_in,
                         const double &abc_tolerance_in,
                         const double &abc_desired_cess_in,
                         const arma::colvec &summary_statistics_scaling_in,
                         const bool &adapt_abc_tolerance_to_cess_in,
                         const bool &adapt_summary_statistics_scaling_in);

  virtual ~ABCLikelihoodEstimator();

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* simulate(const Parameters &parameters) const;

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
};

#endif
