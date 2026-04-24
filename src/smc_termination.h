#ifndef SMCTERMINATION_H
#define SMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"

namespace ilike
{
  /**
   * @file smc_termination.h
   * @brief Defines the SMCTermination class.
   *
   * Implements an SMC termination criterion that terminates based on unconditional. The SMC loop queries this criterion at each iteration to determine whether the algorithm should stop.
   *
   * @namespace ilike
   * @class SMCTermination
   * @brief The smc termination class.
   */


class SMCTermination
{
  
public:
  
  /**
   * @brief Default constructor for SMCTermination.
   */
  SMCTermination();
  /**
   * @brief Destructor for SMCTermination.
   */
  virtual ~SMCTermination();
  
  /**
   * @brief Copy constructor for SMCTermination.
   *
   * @param another The SMCTermination instance to copy from.
   */
  SMCTermination(const SMCTermination &another);
  
  /**
   * @brief Assignment operator for SMCTermination.
   *
   * @param another The SMCTermination instance to copy from.
   */
  void operator=(const SMCTermination &another);
  /**
   * @brief Creates a deep copy of this SMCTermination object.
   *
   * @return The result.
   */
  virtual SMCTermination* duplicate() const=0;
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   *
   * @param score The score.
   *
   * @return The result.
   */
  virtual bool terminate(double score)=0;
  
protected:
  
  /**
   * @brief Copies the state of another SMCTermination into this object.
   *
   * @param another The SMCTermination instance to copy from.
   */
  void make_copy(const SMCTermination &another);
  
};

}

#endif
