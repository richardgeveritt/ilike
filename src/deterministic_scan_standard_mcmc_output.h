#ifndef DETERMINISTICSCANSTANDARDMCMCOUTPUT_H
#define DETERMINISTICSCANSTANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "standard_mcmc_output.h"

class DeterministicScanMCMC;
class MCMC;

class DeterministicScanStandardMCMCOutput : public StandardMCMCOutput
{

public:

  DeterministicScanStandardMCMCOutput();
  
  DeterministicScanStandardMCMCOutput(DeterministicScanMCMC* mcmc_in);

  virtual ~DeterministicScanStandardMCMCOutput();

  DeterministicScanStandardMCMCOutput(const DeterministicScanStandardMCMCOutput &another);

  void operator=(const DeterministicScanStandardMCMCOutput &another);
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
  
  void make_copy(const DeterministicScanStandardMCMCOutput &another);
  
  // Stored here.
  //std::vector<MCMC*> moves;
  //std::vector<size_t> order;
  
  // stored here
  DeterministicScanMCMC* mcmc;

};

#endif
