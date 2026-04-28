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
  /**
   * @file density_likelihood_estimator.h
   * @brief Defines the DensityLikelihoodEstimatorOutput class.
   *
   * Stores the output of a density likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class DensityLikelihoodEstimatorOutput
   * @brief The density likelihood estimator output class.
   */


class DensityLikelihoodEstimatorOutput;
class DensityEstimator;
class DensityLikelihoodEstimatorWorker;
class SequentialDensityLikelihoodEstimatorWorker;
class IndependentProposalKernel;
//class Transform;

class DensityLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the densitylikelihoodestimator operation.
   */
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
  
  /**
   * @brief Performs the ~densitylikelihoodestimator operation.
   */
  virtual ~DensityLikelihoodEstimator();
  
  /**
   * @brief Performs the densitylikelihoodestimator operation.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  DensityLikelihoodEstimator(const DensityLikelihoodEstimator &another);
  
  /**
   * @brief Assignment operator for DensityLikelihoodEstimatorOutput.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const DensityLikelihoodEstimator &another);
  /**
   * @brief Creates a deep copy of this DensityLikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend DensityLikelihoodEstimatorOutput;
  friend DensityLikelihoodEstimatorWorker;
  friend SequentialDensityLikelihoodEstimatorWorker;
  //double evaluate(const Parameters &parameters);
  
  // stored here, to make a deep copy for each Output we make
  /** @brief The density estimator. */
  DensityEstimator* density_estimator;
  /** @brief The subsample density estimator. */
  DensityEstimator* subsample_density_estimator;
  
  // stored here
  /** @brief The proposal. */
  IndependentProposalKernel* proposal;
  /** @brief The subsample proposal. */
  IndependentProposalKernel* subsample_proposal;
  
  //SimulateIndependentProposalPtr simulate_distribution;
  //SimulateIndependentProposalPtr subsample_simulate_distribution;
  
  /** @brief The number of points. */
  size_t number_of_points;
  
  /** @brief The store output. */
  bool store_output;
  
  //Transform* transform;
  
  // Stored here.
  //DensityLikelihoodEstimatorOutput* output;
  
  // Stored here.
  /** @brief The the worker. */
  DensityLikelihoodEstimatorWorker* the_worker;
  
  /** @brief The variables. */
  std::vector<std::string> variables;
  
  /**
   * @brief Copies the state of another DensityLikelihoodEstimatorOutput into this object.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const DensityLikelihoodEstimator &another);
  
  /** @brief HDF5 output file (kept open for the duration of a run). */
  std::shared_ptr<HighFive::File> h5_file;
  /** @brief Path to the HDF5 output file. */
  std::string h5_file_path;
  
};
}

#endif
