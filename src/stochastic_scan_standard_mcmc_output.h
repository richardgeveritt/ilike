#ifndef STOCHASTICSCANSTANDARDMCMCOUTPUT_H
#define STOCHASTICSCANSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
  /**
   * @file stochastic_scan_standard_mcmc_output.h
   * @brief Defines the StochasticScanMCMC class.
   *
   * A stochastic scan MCMC sampler. Implements the MCMC interface using stochastic scan update steps.
   *
   * @namespace ilike
   * @class StochasticScanMCMC
   * @brief The stochastic scan mcmc class.
   */


class StochasticScanMCMC;
class MCMC;

class StochasticScanStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  /**
   * @brief Performs the stochasticscanstandardmcmcoutput operation.
   */
  StochasticScanStandardMCMCOutput();
  
  /**
   * @brief Performs the stochasticscanstandardmcmcoutput operation.
   *
   * @param mcmc_in The mcmc.
   */
  StochasticScanStandardMCMCOutput(StochasticScanMCMC* mcmc_in);
  
  /**
   * @brief Performs the ~stochasticscanstandardmcmcoutput operation.
   */
  virtual ~StochasticScanStandardMCMCOutput();
  
  /**
   * @brief Performs the stochasticscanstandardmcmcoutput operation.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  StochasticScanStandardMCMCOutput(const StochasticScanStandardMCMCOutput &another);
  
  /**
   * @brief Assignment operator for StochasticScanMCMC.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  void operator=(const StochasticScanStandardMCMCOutput &another);
  /**
   * @brief Creates a deep copy of this StochasticScanMCMC object.
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
   * @brief Copies the state of another StochasticScanMCMC into this object.
   *
   * @param another The StochasticScanMCMC instance to copy from.
   */
  void make_copy(const StochasticScanStandardMCMCOutput &another);
  
  // Stored here.
  //std::vector<MCMC*> moves;
  
  //arma::colvec probabilities;
  
  // stored here
  /** @brief The mcmc. */
  StochasticScanMCMC* mcmc;
  
};
}

#endif
