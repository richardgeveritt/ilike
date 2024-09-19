#ifndef METROPOLISHASTINGSSTANDARDMCMCOUTPUT_H
#define METROPOLISHASTINGSSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
class MetropolisHastingsMCMC;

class MetropolisHastingsStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  MetropolisHastingsStandardMCMCOutput();
  
  MetropolisHastingsStandardMCMCOutput(MetropolisHastingsMCMC* mcmc_in);
  
  virtual ~MetropolisHastingsStandardMCMCOutput();
  
  MetropolisHastingsStandardMCMCOutput(const MetropolisHastingsStandardMCMCOutput &another);
  
  void operator=(const MetropolisHastingsStandardMCMCOutput &another);
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
  
  void make_copy(const MetropolisHastingsStandardMCMCOutput &another);
  
  // stored here
  MetropolisHastingsMCMC* mcmc;
  
};
}

#endif
