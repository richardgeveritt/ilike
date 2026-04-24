#ifndef MCMCTERMINATION_H
#define MCMCTERMINATION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file mcmc_termination.h
   * @brief Defines the MCMC class.
   *
   * Provides mcmc functionality.
   *
   * @namespace ilike
   * @class MCMC
   * @brief The mcmc class.
   */


class MCMC;

class MCMCTermination
{
  
public:
  
  /**
   * @brief Performs the mcmctermination operation.
   */
  MCMCTermination();
  /**
   * @brief Performs the ~mcmctermination operation.
   */
  virtual ~MCMCTermination();
  
  /**
   * @brief Performs the mcmctermination operation.
   *
   * @param another The MCMC instance to copy from.
   */
  MCMCTermination(const MCMCTermination &another);
  
  /**
   * @brief Assignment operator for MCMC.
   *
   * @param another The MCMC instance to copy from.
   */
  void operator=(const MCMCTermination &another);
  /**
   * @brief Creates a deep copy of this MCMC object.
   *
   * @return The result.
   */
  virtual MCMCTermination* duplicate() const=0;
  
  /**
   * @brief Returns @c true if the termination criterion has been met.
   *
   * @return The result.
   */
  virtual bool terminate()=0;
  
  /**
   * @brief Sets the parameters.
   *
   * @param mcmc_output The mcmc output.
   */
  virtual void set_parameters(StandardMCMCOutput* mcmc_output)=0;
  
protected:
  
  /**
   * @brief Copies the state of another MCMC into this object.
   *
   * @param another The MCMC instance to copy from.
   */
  void make_copy(const MCMCTermination &another);
  
};
}

#endif
