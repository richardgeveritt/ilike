#ifndef SAMPLEAVERAGEGAUSSIANRECURSIVEPARAMETERESTIMATOR_H
#define SAMPLEAVERAGEGAUSSIANRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "gaussian_recursive_parameter_estimator.h"
#include "particle.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file sample_average_gaussian_recursive_parameter_estimator.h
   * @brief Defines the ProposalKernel class.
   *
   * A generic proposal kernel. Proposes new parameter values during MCMC or SMC moves using a generic distribution centred on the current state.
   *
   * @namespace ilike
   * @class ProposalKernel
   * @brief The proposal kernel class.
   */


class ProposalKernel;

class SampleAverageGaussianRecursiveParameterEstimator : public GaussianRecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Performs the sampleaveragegaussianrecursiveparameterestimator operation.
   */
  SampleAverageGaussianRecursiveParameterEstimator();
  
  /**
   * @brief Performs the ~sampleaveragegaussianrecursiveparameterestimator operation.
   */
  virtual ~SampleAverageGaussianRecursiveParameterEstimator();
  
  /**
   * @brief Performs the sampleaveragegaussianrecursiveparameterestimator operation.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  SampleAverageGaussianRecursiveParameterEstimator(const SampleAverageGaussianRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for ProposalKernel.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void operator=(const SampleAverageGaussianRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this ProposalKernel object.
   *
   * @return The result.
   */
  RecursiveParameterEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  GaussianRecursiveParameterEstimator* gaussian_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  /** @brief The gain. */
  GainPtr gain;
  
  /**
   * @brief Copies the state of another ProposalKernel into this object.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void make_copy(const SampleAverageGaussianRecursiveParameterEstimator &another);
  
};
}

#endif
