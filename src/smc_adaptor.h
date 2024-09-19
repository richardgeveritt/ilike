#ifndef SMCADAPTOR_H
#define SMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
class SMCOutput;
class EnsembleKalmanOutput;

class SMCAdaptor
{
  
public:
  
  SMCAdaptor();
  virtual ~SMCAdaptor();
  
  SMCAdaptor(const SMCAdaptor &another);
  
  void operator=(const SMCAdaptor &another);
  virtual SMCAdaptor* duplicate() const=0;
  
  virtual void smc_adapt(SMCOutput* current_state)=0;
  virtual void ensemble_adapt(EnsembleKalmanOutput* current_state)=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  void make_copy(const SMCAdaptor &another);
  
};
}

#endif
