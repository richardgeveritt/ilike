#ifndef DensityLikelihoodEstimator_H
#define DensityLikelihoodEstimator_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <memory>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

class DensityLikelihoodEstimatorOutput;
class DensityEstimator;
class DensityLikelihoodEstimatorWorker;
class SequentialDensityLikelihoodEstimatorWorker;
class IndependentProposalKernel;

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

  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);

protected:
  
  void specific_change_data(Data* new_data);
  
  friend DensityLikelihoodEstimatorOutput;
  friend DensityLikelihoodEstimatorWorker;
  friend SequentialDensityLikelihoodEstimatorWorker;
  
  // stored here, to make a deep copy for each Output we make
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<DensityEstimator> density_estimator;
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<DensityEstimator> subsample_density_estimator;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<IndependentProposalKernel> proposal;
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<IndependentProposalKernel> subsample_proposal;

  size_t number_of_points;
  
  bool store_output;

  // Stored here.
  //DensityLikelihoodEstimatorOutput* output;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<DensityLikelihoodEstimatorWorker> the_worker;
  
  std::vector<std::string> variables;
  
  void make_copy(const DensityLikelihoodEstimator &another);
  
  std::shared_ptr<HighFive::File> h5_file;
  std::string h5_file_path;

  // values written by DensityLikelihoodEstimatorOutput, flushed per
  // generation: log-likelihoods, times, and the simulated points row by row
  // (point_columns values each)
  std::vector<double> pending_log_likelihoods;
  std::vector<double> pending_times;
  std::vector<double> pending_points;
  size_t point_columns = 0;
  
};
}

#endif
