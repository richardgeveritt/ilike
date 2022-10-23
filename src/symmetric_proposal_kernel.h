#ifndef SYMMETRICPROPOSALKERNEL_H
#define SYMMETRICPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "function_pointers.h"

class SymmetricProposalKernel : public ProposalKernel
{

public:

  SymmetricProposalKernel();
  virtual ~SymmetricProposalKernel();

  SymmetricProposalKernel(const SymmetricProposalKernel &another);

  void operator=(const SymmetricProposalKernel &another);
  //Kernel* duplicate() const;
  virtual SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const=0;
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const SymmetricProposalKernel &another);
  
};

#endif
