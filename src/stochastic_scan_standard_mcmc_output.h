#ifndef STOCHASTICSCANSTANDARDMCMCOUTPUT_H
#define STOCHASTICSCANSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
class StochasticScanMCMC;
class MCMC;

class StochasticScanStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  StochasticScanStandardMCMCOutput();
  
  StochasticScanStandardMCMCOutput(StochasticScanMCMC* mcmc_in);
  
  virtual ~StochasticScanStandardMCMCOutput();
  
  StochasticScanStandardMCMCOutput(const StochasticScanStandardMCMCOutput &another);
  
  void operator=(const StochasticScanStandardMCMCOutput &another);
  MoveOutput* duplicate() const;
  
  /*
   Particle move(RandomNumberGenerator &rng,
   const Particle &particle) const;
   
   Particle subsample_move(RandomNumberGenerator &rng,
   const Particle &particle) const;
   */
  
  MCMC* get_mcmc();
  const MCMC* get_mcmc() const;
  
protected:
  
  //void specific_mcmc_adapt();
  
  void make_copy(const StochasticScanStandardMCMCOutput &another);
  
  // Stored here.
  //std::vector<MCMC*> moves;
  
  //arma::colvec probabilities;
  
  // stored here
  StochasticScanMCMC* mcmc;
  
};
}

#endif
