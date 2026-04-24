#ifndef UNADJUSTEDSTANDARDMCMCOUTPUT_H
#define UNADJUSTEDSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file unadjusted_standard_mcmc_output.h
   * @brief Defines the UnadjustedMCMC class.
   *
   * A unadjusted MCMC sampler. Implements the MCMC interface using unadjusted update steps.
   *
   * @namespace ilike
   * @class UnadjustedMCMC
   * @brief The unadjusted mcmc class.
   */


class UnadjustedMCMC;
class MCMC;

class UnadjustedStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  /**
   * @brief Performs the unadjustedstandardmcmcoutput operation.
   */
  UnadjustedStandardMCMCOutput();
  
  /**
   * @brief Performs the unadjustedstandardmcmcoutput operation.
   *
   * @param mcmc_in The mcmc.
   */
  UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in);
  
  /**
   * @brief Performs the ~unadjustedstandardmcmcoutput operation.
   */
  virtual ~UnadjustedStandardMCMCOutput();
  
  /**
   * @brief Performs the unadjustedstandardmcmcoutput operation.
   *
   * @param another The UnadjustedMCMC instance to copy from.
   */
  UnadjustedStandardMCMCOutput(const UnadjustedStandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for UnadjustedMCMC.
   *
   * @param another The UnadjustedMCMC instance to copy from.
   */
  void operator=(const UnadjustedStandardMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this UnadjustedMCMC object.
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
   * @brief Copies the state of another UnadjustedMCMC into this object.
   *
   * @param another The UnadjustedMCMC instance to copy from.
   */
  void make_copy(const UnadjustedStandardMCMCOutput &another);
  
  // stored here
  /** @brief The mcmc. */
  UnadjustedMCMC* mcmc;
  
};
}

#endif
