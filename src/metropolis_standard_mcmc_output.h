#ifndef METROPOLISSTANDARDMCMCOUTPUT_H
#define METROPOLISSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

class MetropolisMCMC;
class MCMC;

class MetropolisStandardMCMCOutput : public StandardMCMCOutput
{

public:

  MetropolisStandardMCMCOutput();
  
  MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in);

  virtual ~MetropolisStandardMCMCOutput();

  MetropolisStandardMCMCOutput(const MetropolisStandardMCMCOutput &another);

  void operator=(const MetropolisStandardMCMCOutput &another);
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
  
  void make_copy(const MetropolisStandardMCMCOutput &another);
  
  // stored here
  MetropolisMCMC* mcmc;

};

#endif
