#ifndef PCNPROPOSALKERNEL_H
#define PCNPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

// The preconditioned Crank-Nicolson proposal pcn(u, mean, cov, rho) (dev/plan.md
// 3.21, (c1); Cotter, Roberts, Stuart & White 2013): u* = m + rho (u - m) +
// sqrt(1 - rho^2) L xi, with cov = L L^T and xi ~ N(0, I). It leaves the
// Gaussian N(m, cov) invariant (it is reversible with respect to it), so when
// u's prior is that Gaussian the prior cancels from the Metropolis-Hastings
// ratio, which then has only the other factors (e.g. an ABC kernel); its
// acceptance does not degrade with u's dimension. q(u* | u) is
// N(m + rho (u - m), (1 - rho^2) cov): an MH proposal (the prior and q are
// both evaluated, and their ratio cancels numerically).
class PCNProposalKernel : public ProposalKernel
{
public:

  PCNProposalKernel();
  virtual ~PCNProposalKernel();
  PCNProposalKernel(const std::string &variable_in, const arma::colvec &mean_in,
                    const arma::mat &covariance_in, double rho_in);
  PCNProposalKernel(const PCNProposalKernel &another);
  void operator=(const PCNProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  // (C268) simulate, written into the proposal's own storage
  void simulate_into(RandomNumberGenerator &rng,
                     const Particle &particle,
                     ProposalTarget &target) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);

  void make_copy(const PCNProposalKernel &another);
  std::string variable;
  arma::colvec mean;
  double rho;
  arma::mat lower_chol;        // L, cov = L L^T
  arma::mat inverse_step_cov;  // ((1 - rho^2) cov)^-1
  double log_det_step_cov;     // log det((1 - rho^2) cov)
};

}

#endif
