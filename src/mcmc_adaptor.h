#ifndef MCMCADAPTOR_H
#define MCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file mcmc_adaptor.h
   * @brief Defines the MCMCOutput class.
   *
   * Stores and manages the output produced by MCMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class MCMCOutput
   * @brief The mcmc output class.
   */


class MCMCOutput;

class MCMCAdaptor
{
  
public:
  
  /**
   * @brief Performs the mcmcadaptor operation.
   */
  MCMCAdaptor();
  /**
   * @brief Performs the ~mcmcadaptor operation.
   */
  virtual ~MCMCAdaptor();
  
  /**
   * @brief Performs the mcmcadaptor operation.
   *
   * @param another The MCMCOutput instance to copy from.
   */
  MCMCAdaptor(const MCMCAdaptor &another);
  
  /**
   * @brief Assignment operator for MCMCOutput.
   *
   * @param another The MCMCOutput instance to copy from.
   */
  void operator=(const MCMCAdaptor &another);
  /**
   * @brief Creates a deep copy of this MCMCOutput object.
   *
   * @return The result.
   */
  virtual MCMCAdaptor* duplicate() const=0;
  
  void mcmc_adapt(const Particle &latest_particle,
                  size_t iteration_counter);
  
  // MH has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  virtual void specific_mcmc_adapt(const Particle &latest_particle,
                                   size_t iteration_counter)=0;
  
  // not stored here
  /** @brief The proposal. */
  ProposalKernel* proposal;
  
  /**
   * @brief Copies the state of another MCMCOutput into this object.
   *
   * @param another The MCMCOutput instance to copy from.
   */
  void make_copy(const MCMCAdaptor &another);
  
};
}

#endif
