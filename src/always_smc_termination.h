#ifndef ALWAYSSMCTERMINATION_H
#define ALWAYSSMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_termination.h"

namespace ilike
{
  /**
   * @file always_smc_termination.h
   * @brief Defines the AlwaysSMCTermination class.
   *
   * Implements an SMC termination criterion that terminates based on always. The SMC loop queries this criterion at each iteration to determine whether the algorithm should stop.
   *
   * @namespace ilike
   * @class AlwaysSMCTermination
   * @brief An always smc termination derived from SMCTermination.
   */


class AlwaysSMCTermination : public SMCTermination
{
  
public:
  
  /**
   * @brief Default constructor for AlwaysSMCTermination.
   */
  AlwaysSMCTermination();
  
  /**
   * @brief Destructor for AlwaysSMCTermination.
   */
  virtual ~AlwaysSMCTermination();
  
  /**
   * @brief Copy constructor for AlwaysSMCTermination.
   *
   * @param another The AlwaysSMCTermination instance to copy from.
   */
  AlwaysSMCTermination(const AlwaysSMCTermination &another);
  
  /**
   * @brief Assignment operator for AlwaysSMCTermination.
   *
   * @param another The AlwaysSMCTermination instance to copy from.
   */
  void operator=(const AlwaysSMCTermination &another);
  /**
   * @brief Creates a deep copy of this AlwaysSMCTermination object.
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
  
  /**
   * @brief Copies the state of another AlwaysSMCTermination into this object.
   *
   * @param another The AlwaysSMCTermination instance to copy from.
   */
  void make_copy(const AlwaysSMCTermination &another);
  
};
}

#endif
