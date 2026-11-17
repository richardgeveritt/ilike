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
class ProposalKernel;

// Adapts the proposal kernel that owns it from the population of particles,
// between the steps of an SMC algorithm (or of an ensemble Kalman method):
// dev/plan.md 3.2, B2. Called before the particles are resampled and moved,
// on the algorithm's own kernel (which the moves then use read-only).
class SMCAdaptor
{

public:

  SMCAdaptor();
  virtual ~SMCAdaptor();

  SMCAdaptor(const SMCAdaptor &another);
  void operator=(const SMCAdaptor &another);
  virtual SMCAdaptor* duplicate() const=0;

  // The kernel adapted: called when the adaptor is attached, and again on the
  // kernel's copy after the kernel is duplicated.
  virtual void set_proposal(ProposalKernel* proposal_in);

  virtual void smc_adapt(SMCOutput* current_state)=0;
  virtual void ensemble_adapt(EnsembleKalmanOutput* current_state)=0;

protected:

  void make_copy(const SMCAdaptor &another);

  // Not owned: the kernel that owns this adaptor.
  ProposalKernel* proposal;

};

}

#endif
