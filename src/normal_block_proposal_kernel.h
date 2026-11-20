#ifndef NORMALBLOCKPROPOSALKERNEL_H
#define NORMALBLOCKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "hmc_proposal_kernel.h"

namespace ilike
{

// gibbs on a whole vector (stage G2 of the Gibbs design, C216): a draw from
// the block's normal full conditional, as the proposal of a
// Metropolis-Hastings move. When the compile-time check holds (the block's
// own distribution normal, every other term normal with a mean affine in the
// block), log pi(x | rest) = -x' Q x / 2 + b' x + c exactly, so its gradient
// g(x) = b - Q x is linear: from the gradients at the current point x0 and at
// x0 + h e_k (k = 1, ..., d), Q's columns are -(g(x0 + h e_k) - g(x0))/h
// (exact up to rounding, whatever h), and the conditional is
// N(x0 + Q^-1 g(x0), Q^-1). The proposal is drawn from it; q does not depend
// on the current point (given the rest), so the ratio's proposal term is
// log q(x0) - log q(x*), carried by the proposed particle as for HMC (whose
// gradient machinery, index and particle-carried ratio this kernel reuses;
// it has no momentum or step). The ratio is 1 up to rounding and is
// computed, so the chain is exact whatever the rounding (and, if the target
// were not Gaussian, this would still be a valid independence proposal).
class NormalBlockProposalKernel : public HMCProposalKernel
{

public:

  NormalBlockProposalKernel();
  virtual ~NormalBlockProposalKernel();
  NormalBlockProposalKernel(const std::string &variable_in);

  NormalBlockProposalKernel(const NormalBlockProposalKernel &another);
  void operator=(const NormalBlockProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  Parameters trajectory(RandomNumberGenerator &rng,
                        const Particle &particle,
                        double &log_ratio) const;

  // the conditional's mean and the Cholesky factor U of its precision
  // (Q = U' U), from the particle's point; false if Q is not positive definite
  bool conditional(const Particle &particle,
                   const Parameters &current,
                   arma::colvec &mean,
                   arma::mat &upper) const;
};

}

#endif
