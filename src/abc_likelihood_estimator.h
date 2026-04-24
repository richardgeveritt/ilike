#ifndef ABCLIKELIHOODESTIMATOR_H
#define ABCLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "parameters.h"
#include "abc_kernel_factor.h"

namespace ilike
{
  /**
   * @file abc_likelihood_estimator.h
   * @brief Defines the ABCLikelihoodEstimatorOutput class.
   *
   * Provides functions for an abc likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class ABCLikelihoodEstimatorOutput
   * @brief The abc likelihood estimator output class.
   */



class ABCLikelihoodEstimatorOutput;
class IndependentProposalKernel;

class ABCLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the abclikelihoodestimator operation.
   */
  ABCLikelihoodEstimator();
  
  ABCLikelihoodEstimator(RandomNumberGenerator* rng_in,
                         size_t* seed_in,
                         Data* data_in,
                         ABCKernelFactor* abc_kernel_in,
                         bool smcfixed_flag_in);
  
  /**
   * @brief Performs the ~abclikelihoodestimator operation.
   */
  virtual ~ABCLikelihoodEstimator();
  
  /**
   * @brief Performs the abclikelihoodestimator operation.
   *
   * @param another The ABCLikelihoodEstimatorOutput instance to copy from.
   */
  ABCLikelihoodEstimator(const ABCLikelihoodEstimator &another);
  
  /**
   * @brief Assignment operator for ABCLikelihoodEstimatorOutput.
   *
   * @param another The ABCLikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const ABCLikelihoodEstimator &another);
  /**
   * @brief Creates a deep copy of this ABCLikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
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
  
private:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend ABCLikelihoodEstimatorOutput;
  
  // Stored here.
  /** @brief The abc kernel. */
  ABCKernelFactor* abc_kernel;
  
  /** @brief The log likelihood file stream. */
  std::ofstream log_likelihood_file_stream;
  
  /**
   * @brief Copies the state of another ABCLikelihoodEstimatorOutput into this object.
   *
   * @param another The ABCLikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const ABCLikelihoodEstimator &another);
  
};
}

#endif
