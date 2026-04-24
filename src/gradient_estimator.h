#ifndef GRADIENTESTIMATOR_H
#define GRADIENTESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include "particle.h"
#include "data_subsampler.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file gradient_estimator.h
   * @brief Defines the ReinforceGradientEstimator class.
   *
   * Estimates the reinforce gradient for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class ReinforceGradientEstimator
   * @brief The reinforce gradient estimator class.
   */



class ReinforceGradientEstimator;
class GradientEstimatorOutput;
class ProposalKernel;
class DirectGradientEstimatorOutput;

class GradientEstimator
{
  
public:
  
  /**
   * @brief Performs the gradientestimator operation.
   */
  GradientEstimator();
  /**
   * @brief Performs the ~gradientestimator operation.
   */
  virtual ~GradientEstimator();
  
  /**
   * @brief Performs the gradientestimator operation.
   *
   * @param proposal_in The proposal.
   */
  GradientEstimator(const ProposalKernel* proposal_in);
  
  /**
   * @brief Performs the gradientestimator operation.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  GradientEstimator(const GradientEstimator &another);
  
  /**
   * @brief Assignment operator for ReinforceGradientEstimator.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void operator=(const GradientEstimator &another);
  /**
   * @brief Creates a deep copy of this ReinforceGradientEstimator object.
   *
   * @return The result.
   */
  virtual GradientEstimator* duplicate() const=0;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  virtual GradientEstimatorOutput* initialise()=0;
  
  /**
   * @brief Sets the proposal.
   *
   * @param proposal_in The proposal.
   */
  void set_proposal(const ProposalKernel* proposal_in);
  
  /*
   virtual arma::mat get_gradient_of_log(const std::string &variable,
   Particle &particle)=0;
   
   virtual arma::mat get_gradient_of_log(const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters)=0;
   
   virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
   Particle &particle)=0;
   
   virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters)=0;
   */
  
protected:
  
  friend ReinforceGradientEstimator;
  friend DirectGradientEstimatorOutput;
  // not stored here. A pointer to the proposal in which this gradient estimator is contained
  /** @brief The proposal. */
  const ProposalKernel* proposal;
  
  // Not stored here. Stored in "main'.
  /** @brief The rng. */
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  /** @brief The seed. */
  size_t* seed;
  
  /**
   * @brief Copies the state of another ReinforceGradientEstimator into this object.
   *
   * @param another The ReinforceGradientEstimator instance to copy from.
   */
  void make_copy(const GradientEstimator &another);
  
};
}

#endif
