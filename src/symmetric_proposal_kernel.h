#ifndef SYMMETRICPROPOSALKERNEL_H
#define SYMMETRICPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file symmetric_proposal_kernel.h
   * @brief Defines the SymmetricProposalKernel class.
   *
   * A symmetric proposal kernel. Proposes new parameter values during MCMC or SMC moves using a symmetric distribution centred on the current state.
   *
   * @namespace ilike
   * @class SymmetricProposalKernel
   * @brief A symmetric proposal kernel derived from ProposalKernel.
   */


class SymmetricProposalKernel : public ProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for SymmetricProposalKernel.
   */
  SymmetricProposalKernel();
  /**
   * @brief Destructor for SymmetricProposalKernel.
   */
  virtual ~SymmetricProposalKernel();
  
  /**
   * @brief Copy constructor for SymmetricProposalKernel.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  SymmetricProposalKernel(const SymmetricProposalKernel &another);
  
  /**
   * @brief Assignment operator for SymmetricProposalKernel.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  void operator=(const SymmetricProposalKernel &another);
  //Kernel* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a symmetric_proposal_kernel pointer.
   *
   * @return The result.
   */
  virtual SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  /**
   * @brief Copies the state of another SymmetricProposalKernel into this object.
   *
   * @param another The SymmetricProposalKernel instance to copy from.
   */
  void make_copy(const SymmetricProposalKernel &another);
  
};
}

#endif
