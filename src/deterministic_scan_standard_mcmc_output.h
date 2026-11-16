#ifndef DETERMINISTICSCANSTANDARDMCMCOUTPUT_H
#define DETERMINISTICSCANSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"
#include "optional_delete.h"

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

  // owns_mcmc_in false: share the kernel (not deleted or duplicated here;
  // it must outlive this output, as the algorithm that owns it does).
  DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in, bool owns_mcmc_in);
  
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
  // The move. Owned (deleted with this output, duplicated with it) for a move
  // that adapts: each chain adapts its own copy. Not owned when the move is
  // shared between the particles' outputs (dev/plan.md, owner's decision 2,
  // 2026-11-15); then it must outlive this output, as the algorithm that owns
  // it does. Whether it is owned is mcmc.get_deleter().owns.
  std::unique_ptr<DeterministicScanMCMC, OptionalDelete> mcmc;
  
};
}

#endif
