#ifndef REINFORCEGRADIENTESTIMATOR_H
#define REINFORCEGRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gradient_estimator.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
  /**
   * @file reinforce_gradient_estimator.h
   * @brief Defines the GradientEstimatorOutput class.
   *
   * Stores the output of a gradient estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class GradientEstimatorOutput
   * @brief The gradient estimator output class.
   */


class GradientEstimatorOutput;

class ReinforceGradientEstimator : public GradientEstimator
{
  
public:
  
  /**
   * @brief Performs the reinforcegradientestimator operation.
   */
  ReinforceGradientEstimator();
  
  /**
   * @brief Performs the reinforcegradientestimator operation.
   *
   * @param subsampler_in The subsampler.
   */
  ReinforceGradientEstimator(DataSubsampler* subsampler_in);
  
  /**
   * @brief Performs the ~reinforcegradientestimator operation.
   */
  virtual ~ReinforceGradientEstimator();
  
  /**
   * @brief Performs the reinforcegradientestimator operation.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  ReinforceGradientEstimator(const ReinforceGradientEstimator &another);
  
  /**
   * @brief Assignment operator for GradientEstimatorOutput.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void operator=(const ReinforceGradientEstimator &another);
  /**
   * @brief Creates a deep copy of this GradientEstimatorOutput object.
   *
   * @return The result.
   */
  GradientEstimator* duplicate() const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  GradientEstimatorOutput* initialise();
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const std::vector<arma::mat> &auxiliary_variables,
                                const Index* index,
                                const Particle &particle);
  
  /*
   arma::mat get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const std::vector<arma::mat> &auxiliary_variables,
                                          const Index* index,
                                          const Particle &particle);
  
  /**
   * @brief Simulates auxiliary variables.
   *
   * @return The result.
   */
  boost::unordered_map< std::string, std::vector<arma::mat>> simulate_auxiliary_variables();
  
  /*
   arma::mat subsample_get_gradient_of_log(const std::string &variable,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
protected:
  
  /**
   * @brief Copies the state of another GradientEstimatorOutput into this object.
   *
   * @param another The GradientEstimatorOutput instance to copy from.
   */
  void make_copy(const ReinforceGradientEstimator &another);
  
  // can be generalised, but let's not worry about that for now
  /** @brief The gaussian proposal. */
  GaussianIndependentProposalKernel gaussian_proposal;
  
  /** @brief The num points. */
  size_t num_points;
  /** @brief The size of subsample. */
  size_t size_of_subsample;
  
  // Not stored here. Stored in "main'.
  /** @brief The subsampler. */
  DataSubsampler* subsampler;
  
};
}

#endif
