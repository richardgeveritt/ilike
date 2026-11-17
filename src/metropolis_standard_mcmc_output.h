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

class MetropolisMCMC;
class MCMC;

class MetropolisStandardMCMCOutput : public StandardMCMCOutput
{
  
public:
  
  MetropolisStandardMCMCOutput();
  
  MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in);

  // owns_mcmc_in false: share the kernel (not deleted or duplicated here;
  // it must outlive this output, as the algorithm that owns it does).
  MetropolisStandardMCMCOutput(MetropolisMCMC* mcmc_in, bool owns_mcmc_in);
  
  virtual ~MetropolisStandardMCMCOutput();
  
  MetropolisStandardMCMCOutput(const MetropolisStandardMCMCOutput &another);
  
  void operator=(const MetropolisStandardMCMCOutput &another);
  MoveOutput* duplicate() const;

  MCMC* get_mcmc();
  const MCMC* get_mcmc() const;
  
protected:

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
