#ifndef CONJUGATEBLOCKPROPOSALKERNEL_H
#define CONJUGATEBLOCKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "hmc_proposal_kernel.h"

namespace ilike
{

// gibbs on a Dirichlet-distributed probability vector p with categorical(p)
// and multinom(n, p) terms (Gibbs design item 5, C228): its log full
// conditional is exactly sum_k c_k log p_k + const on the simplex. From the
// log density at p and at the K points p^(j) = p e^(h e_j)/(1 + p_j (e^h - 1))
// (p_j scaled by e^h, renormalised): with D_j the differences from p and
// L_j = log(1 + p_j (e^h - 1)), D_j = c_j h - S L_j, S = sum c_k, so
// S = sum D/(h - sum L) and c_j = (D_j + S L_j)/h; the conditional is
// Dirichlet(c + 1). Drawn as the proposal of a Metropolis-Hastings move (q
// independent of the current p given the rest: the particle carries
// log q(p) - log q(p*), as HMC's kernel energy; ratio 1 up to rounding).
class DirichletBlockProposalKernel : public HMCProposalKernel
{
public:
  DirichletBlockProposalKernel();
  virtual ~DirichletBlockProposalKernel();
  DirichletBlockProposalKernel(const std::string &variable_in);
  DirichletBlockProposalKernel(const DirichletBlockProposalKernel &another);
  void operator=(const DirichletBlockProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  Parameters trajectory(RandomNumberGenerator &rng, const Particle &particle, double &log_ratio) const;
};

// gibbs([mu, v]) on a mean and a variance with a normal-inverse-gamma full
// conditional (C228): v ~ invgamma, mu ~ norm(m0, c sqrt(v)), and every
// other term norm(a + b mu, c' sqrt(v)). The log density is exactly
// alpha log v + (beta + gamma mu + delta mu^2)/v + const: fitted from the log
// density at five points around the current one (a 5 x 5 solve, exact up
// to rounding), then v drawn from its marginal, inverse gamma with shape
// -alpha - 3/2 and scale -(beta - gamma^2/(4 delta)), and mu given v from
// N(-gamma/(2 delta), -v/(2 delta)); as a Metropolis-Hastings proposal, as
// above. mean_first: whether mu is the first variable of the block.
class NormalInverseGammaBlockProposalKernel : public HMCProposalKernel
{
public:
  NormalInverseGammaBlockProposalKernel();
  virtual ~NormalInverseGammaBlockProposalKernel();
  NormalInverseGammaBlockProposalKernel(const std::string &variable_in, bool mean_first_in);
  NormalInverseGammaBlockProposalKernel(const NormalInverseGammaBlockProposalKernel &another);
  void operator=(const NormalInverseGammaBlockProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  Parameters trajectory(RandomNumberGenerator &rng, const Particle &particle, double &log_ratio) const;
protected:
  bool mean_first;
};

}

#endif
