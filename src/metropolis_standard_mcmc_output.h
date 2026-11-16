#ifndef METROPOLISSTANDARDMCMCOUTPUT_H
#define METROPOLISSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"
#include "optional_delete.h"

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

  // owns_mcmc_in false: share the kernel (not deleted or duplicated here;
  // it must outlive this output, as the algorithm that owns it does).
  MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in, bool owns_mcmc_in);
  
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
  // The move. Owned (deleted with this output, duplicated with it) for a move
  // that adapts: each chain adapts its own copy. Not owned when the move is
  // shared between the particles' outputs (dev/plan.md, owner's decision 2,
  // 2026-11-15); then it must outlive this output, as the algorithm that owns
  // it does. Whether it is owned is mcmc.get_deleter().owns.
  std::unique_ptr<MetropolisMCMC, OptionalDelete> mcmc;
  
};
}

#endif
