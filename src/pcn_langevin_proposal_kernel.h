#ifndef PCNLANGEVINPROPOSALKERNEL_H
#define PCNLANGEVINPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "langevin_proposal_kernel.h"

namespace ilike
{

// pCN-MALA (the preconditioned Crank-Nicolson Langevin proposal, Cotter,
// Roberts, Stuart & White 2013; dev/plan.md 3.21 (c3)): for u with a
// Gaussian distribution N(m, C) and the other factors exp(-Phi(u)) (e.g. an
// ABC kernel with a gradient),
//   u* = m + a (u - m) - b C DPhi(u) + c L xi,  xi ~ N(0, I), C = L L^T,
// a = (2 - delta)/(2 + delta), b = 2 delta/(2 + delta), c = sqrt(8 delta)/(2 + delta).
// DPhi(u) = -grad log pi(u) - C^-1 (u - m) (the gradient of the move's whole
// target, less that of the Gaussian). With DPhi = 0 it is pCN (rho = a). An MH
// proposal: q(u* | u) = N(u*; mean, c^2 C) is in the ratio, so it is exact
// for any delta.
class PCNLangevinProposalKernel : public LangevinProposalKernel
{
public:

  PCNLangevinProposalKernel();
  virtual ~PCNLangevinProposalKernel();
  PCNLangevinProposalKernel(const std::string &variable_in, const arma::colvec &mean_in,
                            const arma::mat &covariance_in, double delta_in,
                            GradientEstimator* gradient_estimator_in);
  PCNLangevinProposalKernel(const PCNLangevinProposalKernel &another);
  void operator=(const PCNLangevinProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  void simulate_into(RandomNumberGenerator &rng,
                     const Particle &particle,
                     ProposalTarget &target) const;

  // the proposal's mean at the particle
  arma::colvec proposal_mean(const Particle &particle) const;

  void make_copy(const PCNLangevinProposalKernel &another);

  std::string variable;
  arma::colvec mean;
  arma::mat covariance;
  arma::mat inverse_covariance;
  arma::mat lower_chol;
  double a = 1.0, b = 0.0, c = 0.0;
  arma::mat inverse_step_covariance; // (c^2 C)^-1
  double log_det_step_covariance = 0.0;
};

}

#endif
