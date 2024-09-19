#ifndef MCMCADAPTOR_H
#define MCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
class MCMCOutput;

class MCMCAdaptor
{
  
public:
  
  MCMCAdaptor();
  virtual ~MCMCAdaptor();
  
  MCMCAdaptor(const MCMCAdaptor &another);
  
  void operator=(const MCMCAdaptor &another);
  virtual MCMCAdaptor* duplicate() const=0;
  
  void mcmc_adapt(const Particle &latest_particle,
                  size_t iteration_counter);
  
  // MH has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  virtual void specific_mcmc_adapt(const Particle &latest_particle,
                                   size_t iteration_counter)=0;
  
  // not stored here
  ProposalKernel* proposal;
  
  void make_copy(const MCMCAdaptor &another);
  
};
}

#endif
