#ifndef UNADJUSTEDSTANDARDMCMCOUTPUT_H
#define UNADJUSTEDSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

namespace ilike
{
class UnadjustedMCMC;
class MCMC;

class UnadjustedStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  UnadjustedStandardMCMCOutput();
  
  UnadjustedStandardMCMCOutput(UnadjustedMCMC* mcmc_in);
  
  virtual ~UnadjustedStandardMCMCOutput();
  
  UnadjustedStandardMCMCOutput(const UnadjustedStandardMCMCOutput &another);
  
  void operator=(const UnadjustedStandardMCMCOutput &another);
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
  
  void make_copy(const UnadjustedStandardMCMCOutput &another);
  
  // stored here
  UnadjustedMCMC* mcmc;
  
};
}

#endif
