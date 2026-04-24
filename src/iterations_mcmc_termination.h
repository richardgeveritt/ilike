#ifndef ITERATIONSMCMCTERMINATION_H
#define ITERATIONSMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

// Checks to see

namespace ilike
{
  /**
   * @file iterations_mcmc_termination.h
   * @brief Defines the IterationsMCMCTermination class.
   *
   * Implements an MCMC termination criterion based on iterations. The sampler queries this criterion after each step.
   *
   * @namespace ilike
   * @class IterationsMCMCTermination
   * @brief An iterations mcmc termination derived from MCMCTermination.
   */


class IterationsMCMCTermination : public MCMCTermination
{
  
public:
  
  /**
   * @brief Default constructor for IterationsMCMCTermination.
   */
  IterationsMCMCTermination();
  
  /**
   * @brief Constructs a IterationsMCMCTermination object.
   *
   * @param number_of_iterations_in The number of iterations.
   */
  IterationsMCMCTermination(size_t number_of_iterations_in);
  
  /**
   * @brief Destructor for IterationsMCMCTermination.
   */
  virtual ~IterationsMCMCTermination();
  
  /**
   * @brief Copy constructor for IterationsMCMCTermination.
   *
   * @param another The IterationsMCMCTermination instance to copy from.
   */
  IterationsMCMCTermination(const IterationsMCMCTermination &another);
  
  /**
   * @brief Assignment operator for IterationsMCMCTermination.
   *
   * @param another The IterationsMCMCTermination instance to copy from.
   */
  void operator=(const IterationsMCMCTermination &another);
  /**
   * @brief Creates a deep copy of this IterationsMCMCTermination object.
   *
   * @return The result.
   */
  MCMCTermination* duplicate() const;
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   *
   * @return The result.
   */
  bool terminate();
  
  /**
   * @brief Sets the parameters.
   *
   * @param mcmc_output The mcmc output.
   */
  void set_parameters(StandardMCMCOutput* mcmc_output);
  
protected:
  
  // not stored here
  /** @brief The counter. */
  size_t* counter;
  
  /** @brief The number of iterations. */
  size_t number_of_iterations;
  
  /**
   * @brief Copies the state of another IterationsMCMCTermination into this object.
   *
   * @param another The IterationsMCMCTermination instance to copy from.
   */
  void make_copy(const IterationsMCMCTermination &another);
  
};
}

#endif
