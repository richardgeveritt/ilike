#ifndef DensityLikelihoodEstimator_H
#define DensityLikelihoodEstimator_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"

class DensityLikelihoodEstimatorOutput;
class DensityEstimator;
class DensityLikelihoodEstimatorWorker;
class SequentialDensityLikelihoodEstimatorWorker;

class DensityLikelihoodEstimator : public LikelihoodEstimator
{

public:

  DensityLikelihoodEstimator();

  // Takes control of DensityEstimator memory.
  DensityLikelihoodEstimator(RandomNumberGenerator* rng_in,
                             size_t* seed_in,
                             Data* data_in,
                             size_t number_of_points_in,
                             bool smcfixed_flag_in,
                             DensityEstimator* density_estimator_in,
                             SimulateIndependentProposalPtr simulate_distribution_in,
                             bool parallel_in);
  
  virtual ~DensityLikelihoodEstimator();

  DensityLikelihoodEstimator(const DensityLikelihoodEstimator &another);

  void operator=(const DensityLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

protected:

  friend DensityLikelihoodEstimatorOutput;
  friend DensityLikelihoodEstimatorWorker;
  friend SequentialDensityLikelihoodEstimatorWorker;
  //double evaluate(const Parameters &parameters);
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;

  // Stored here.
  DensityEstimator* density_estimator;
  DensityEstimator* subsample_density_estimator;
  
  SimulateIndependentProposalPtr simulate_distribution;
  SimulateIndependentProposalPtr subsample_simulate_distribution;
  
  size_t number_of_points;

  // Stored here.
  //DensityLikelihoodEstimatorOutput* output;
  
  // Stored here.
  DensityLikelihoodEstimatorWorker* the_worker;

  void make_copy(const DensityLikelihoodEstimator &another);

};

#endif
