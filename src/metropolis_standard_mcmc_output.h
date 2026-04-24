#ifndef METROPOLISSTANDARDMCMCOUTPUT_H
#define METROPOLISSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file metropolis_standard_mcmc_output.h
   * @brief Defines the MetropolisMCMC class.
   *
   * A metropolis MCMC sampler. Implements the MCMC interface using metropolis update steps.
   *
   * @namespace ilike
   * @class MetropolisMCMC
   * @brief The metropolis mcmc class.
   */


class MetropolisMCMC;
class MCMC;

class MetropolisStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  /**
   * @brief Performs the metropolisstandardmcmcoutput operation.
   */
  MetropolisStandardMCMCOutput();
  
  /**
   * @brief Performs the metropolisstandardmcmcoutput operation.
   *
   * @param mcmc_in The mcmc.
   */
  MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in);
  
  /**
   * @brief Performs the ~metropolisstandardmcmcoutput operation.
   */
  virtual ~MetropolisStandardMCMCOutput();
  
  /**
   * @brief Performs the metropolisstandardmcmcoutput operation.
   *
   * @param another The MetropolisMCMC instance to copy from.
   */
  MetropolisStandardMCMCOutput(const MetropolisStandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for MetropolisMCMC.
   *
   * @param another The MetropolisMCMC instance to copy from.
   */
  void operator=(const MetropolisStandardMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this MetropolisMCMC object.
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
   * @brief Copies the state of another MetropolisMCMC into this object.
   *
   * @param another The MetropolisMCMC instance to copy from.
   */
  void make_copy(const MetropolisStandardMCMCOutput &another);
  
  // stored here
  /** @brief The mcmc. */
  MetropolisMCMC* mcmc;
  
};
}

#endif
