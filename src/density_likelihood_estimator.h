#ifndef DensityLikelihoodEstimator_H
#define DensityLikelihoodEstimator_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"

class DensityLikelihoodEstimatorOutput;
class DensityEstimator;
class DensityLikelihoodEstimatorWorker;
class SequentialDensityLikelihoodEstimatorWorker;
class IndependentProposalKernel;
//class Transform;

class DensityLikelihoodEstimator : public LikelihoodEstimator
{

public:

  DensityLikelihoodEstimator();

  DensityLikelihoodEstimator(RandomNumberGenerator* rng_in,
                             size_t* seed_in,
                             Data* data_in,
                             const Parameters &algorithm_parameters_in,
                             size_t number_of_points_in,
                             bool smcfixed_flag_in,
                             DensityEstimator* density_estimator_in,
                             IndependentProposalKernel* proposal_in,
                             bool make_subsample_version_in,
                             bool store_output_in,
                             bool parallel_in,
                             size_t grain_size_in);
  
  virtual ~DensityLikelihoodEstimator();

  DensityLikelihoodEstimator(const DensityLikelihoodEstimator &another);

  void operator=(const DensityLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);

protected:
  
  void specific_change_data(Data* new_data);

  friend DensityLikelihoodEstimatorOutput;
  friend DensityLikelihoodEstimatorWorker;
  friend SequentialDensityLikelihoodEstimatorWorker;
  //double evaluate(const Parameters &parameters);
  
  // stored here, to make a deep copy for each Output we make
  DensityEstimator* density_estimator;
  DensityEstimator* subsample_density_estimator;
  
  // stored here
  IndependentProposalKernel* proposal;
  IndependentProposalKernel* subsample_proposal;
  
  //SimulateIndependentProposalPtr simulate_distribution;
  //SimulateIndependentProposalPtr subsample_simulate_distribution;
  
  size_t number_of_points;
  
  bool store_output;

  //Transform* transform;
  
  // Stored here.
  //DensityLikelihoodEstimatorOutput* output;
  
  // Stored here.
  DensityLikelihoodEstimatorWorker* the_worker;
  
  std::vector<std::string> variables;

  void make_copy(const DensityLikelihoodEstimator &another);
  
  std::ofstream log_likelihood_file_stream;
  std::ofstream file_stream;
  std::ofstream time_file_stream;

};

#endif
