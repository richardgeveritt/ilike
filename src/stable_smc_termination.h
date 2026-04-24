#ifndef STABLESMCTERMINATION_H
#define STABLESMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_termination.h"

// Checks to see

namespace ilike
{
  /**
   * @file stable_smc_termination.h
   * @brief Defines the StableSMCTermination class.
   *
   * Implements an SMC termination criterion that terminates based on stable. The SMC loop queries this criterion at each iteration to determine whether the algorithm should stop.
   *
   * @namespace ilike
   * @class StableSMCTermination
   * @brief A stable smc termination derived from SMCTermination.
   */


class StableSMCTermination : public SMCTermination
{
  
public:
  
  /**
   * @brief Default constructor for StableSMCTermination.
   */
  StableSMCTermination();
  
  StableSMCTermination(size_t number_in_a_row_in,
                       double threshold_in);
  
  /**
   * @brief Destructor for StableSMCTermination.
   */
  virtual ~StableSMCTermination();
  
  /**
   * @brief Copy constructor for StableSMCTermination.
   *
   * @param another The StableSMCTermination instance to copy from.
   */
  StableSMCTermination(const StableSMCTermination &another);
  
  /**
   * @brief Assignment operator for StableSMCTermination.
   *
   * @param another The StableSMCTermination instance to copy from.
   */
  void operator=(const StableSMCTermination &another);
  /**
   * @brief Creates a deep copy of this StableSMCTermination object.
   *
   * @return The result.
   */
  SMCTermination* duplicate() const;
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   *
   * @param score The score.
   *
   * @return The result.
   */
  bool terminate(double score);
  
protected:
  
  /** @brief The counter. */
  size_t counter;
  
  /** @brief The number in a row. */
  size_t number_in_a_row;
  /** @brief The threshold. */
  double threshold;
  
  /**
   * @brief Copies the state of another StableSMCTermination into this object.
   *
   * @param another The StableSMCTermination instance to copy from.
   */
  void make_copy(const StableSMCTermination &another);
  
};
}

#endif
