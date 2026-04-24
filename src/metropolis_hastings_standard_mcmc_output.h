#ifndef METROPOLISHASTINGSSTANDARDMCMCOUTPUT_H
#define METROPOLISHASTINGSSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file metropolis_hastings_standard_mcmc_output.h
   * @brief Defines the MetropolisHastingsMCMC class.
   *
   * A metropolis hastings MCMC sampler. Implements the MCMC interface using metropolis hastings update steps.
   *
   * @namespace ilike
   * @class MetropolisHastingsMCMC
   * @brief The metropolis hastings mcmc class.
   */


class MetropolisHastingsMCMC;

class MetropolisHastingsStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  /**
   * @brief Performs the metropolishastingsstandardmcmcoutput operation.
   */
  MetropolisHastingsStandardMCMCOutput();
  
  /**
   * @brief Performs the metropolishastingsstandardmcmcoutput operation.
   *
   * @param mcmc_in The mcmc.
   */
  MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in);
  
  /**
   * @brief Performs the ~metropolishastingsstandardmcmcoutput operation.
   */
  virtual ~MetropolisHastingsStandardMCMCOutput();
  
  /**
   * @brief Performs the metropolishastingsstandardmcmcoutput operation.
   *
   * @param another The MetropolisHastingsMCMC instance to copy from.
   */
  MetropolisHastingsStandardMCMCOutput(const MetropolisHastingsStandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for MetropolisHastingsMCMC.
   *
   * @param another The MetropolisHastingsMCMC instance to copy from.
   */
  void operator=(const MetropolisHastingsStandardMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this MetropolisHastingsMCMC object.
   *
   * @return The result.
   */
  MoveOutput* duplicate() const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   const Particle &particle) const;
   
   Particle subsample_move(RandomNumberGenerator &rng,
   const Particle &particle) const;
   */
  
  MCMC* get_mcmc();
  /**
   * @brief Returns the mcmc.
   *
   * @return The result.
   */
  const MCMC* get_mcmc() const;
  
protected:
  
  //void specific_mcmc_adapt();
  
  /**
   * @brief Copies the state of another MetropolisHastingsMCMC into this object.
   *
   * @param another The MetropolisHastingsMCMC instance to copy from.
   */
  void make_copy(const MetropolisHastingsStandardMCMCOutput &another);
  
  // stored here
  /** @brief The mcmc. */
  MetropolisHastingsMCMC* mcmc;
  
};
}

#endif
