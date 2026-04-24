#ifndef SEMCMCTERMINATION_H
#define SEMCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "mcmc_termination.h"

namespace ilike
{
  /**
   * @file se_mcmc_termination.h
   * @brief Defines the SEMCMCTermination class.
   *
   * Implements an MCMC termination criterion based on se. The sampler queries this criterion after each step.
   *
   * @namespace ilike
   * @class SEMCMCTermination
   * @brief A semcmc termination derived from MCMCTermination.
   */


class SEMCMCTermination : public MCMCTermination
{
  
public:
  
  /**
   * @brief Default constructor for SEMCMCTermination.
   */
  SEMCMCTermination();
  
  SEMCMCTermination(double threshold_in,
                    size_t max_number_of_iterations_in);
  
  /**
   * @brief Destructor for SEMCMCTermination.
   */
  virtual ~SEMCMCTermination();
  
  /**
   * @brief Copy constructor for SEMCMCTermination.
   *
   * @param another The SEMCMCTermination instance to copy from.
   */
  SEMCMCTermination(const SEMCMCTermination &another);
  
  /**
   * @brief Assignment operator for SEMCMCTermination.
   *
   * @param another The SEMCMCTermination instance to copy from.
   */
  void operator=(const SEMCMCTermination &another);
  /**
   * @brief Creates a deep copy of this SEMCMCTermination object.
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
  
  /** @brief The threshold. */
  double threshold;
  
  /** @brief The max number of iterations. */
  size_t max_number_of_iterations;
  
  /**
   * @brief Copies the state of another SEMCMCTermination into this object.
   *
   * @param another The SEMCMCTermination instance to copy from.
   */
  void make_copy(const SEMCMCTermination &another);
  
};
}

#endif
