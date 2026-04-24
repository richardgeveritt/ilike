#ifndef DETERMINISTICSCANSTANDARDMCMCOUTPUT_H
#define DETERMINISTICSCANSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file deterministic_scan_standard_mcmc_output.h
   * @brief Defines the DeterministicScanMCMC class.
   *
   * A deterministic scan MCMC sampler. Implements the MCMC interface using deterministic scan update steps.
   *
   * @namespace ilike
   * @class DeterministicScanMCMC
   * @brief The deterministic scan mcmc class.
   */


class DeterministicScanMCMC;
class MCMC;

class DeterministicScanStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  /**
   * @brief Performs the deterministicscanstandardmcmcoutput operation.
   */
  DeterministicScanStandardMCMCOutput();
  
  /**
   * @brief Performs the deterministicscanstandardmcmcoutput operation.
   *
   * @param mcmc_in The mcmc.
   */
  DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in);
  
  /**
   * @brief Performs the ~deterministicscanstandardmcmcoutput operation.
   */
  virtual ~DeterministicScanStandardMCMCOutput();
  
  /**
   * @brief Performs the deterministicscanstandardmcmcoutput operation.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  DeterministicScanStandardMCMCOutput(const DeterministicScanStandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for DeterministicScanMCMC.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  void operator=(const DeterministicScanStandardMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this DeterministicScanMCMC object.
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
   * @brief Copies the state of another DeterministicScanMCMC into this object.
   *
   * @param another The DeterministicScanMCMC instance to copy from.
   */
  void make_copy(const DeterministicScanStandardMCMCOutput &another);
  
  // Stored here.
  //std::vector<MCMC*> moves;
  //std::vector<size_t> order;
  
  // stored here
  /** @brief The mcmc. */
  DeterministicScanMCMC* mcmc;
  
};
}

#endif
