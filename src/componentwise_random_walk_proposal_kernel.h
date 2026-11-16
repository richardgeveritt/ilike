#ifndef COMPONENTWISERANDOMWALKPROPOSALKERNEL_H
#define COMPONENTWISERANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "gaussian_random_walk_proposal_kernel.h"

namespace ilike
{

// Random-scan componentwise Gaussian random walk (Andrieu & Thoms 2008,
// Algorithm 5; dev/plan.md 3.2, A3): one component k of the variable, chosen
// uniformly, moves by N(0, lambda_k [Sigma]_kk); the others stay. Sigma is the
// covariance stored as in GaussianRandomWalkProposalKernel (only its diagonal
// is used here); lambda_k are per-component scales, all starting at the
// kernel's scale. The move is symmetric, so in a Metropolis move its density
// cancels and is never needed (it is singular: refused if asked for).
class ComponentwiseRandomWalkProposalKernel : public GaussianRandomWalkProposalKernel
{

public:

  ComponentwiseRandomWalkProposalKernel();
  virtual ~ComponentwiseRandomWalkProposalKernel();

  ComponentwiseRandomWalkProposalKernel(const std::string &variable_name_in,
                                        const arma::mat &covariance_in,
                                        double scale_in);

  ComponentwiseRandomWalkProposalKernel(const ComponentwiseRandomWalkProposalKernel &another);
  void operator=(const ComponentwiseRandomWalkProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;

  // lambda_1, ..., lambda_d, for the componentwise adaptor
  arma::colvec* get_component_scales(const std::string &variable);

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                  const Parameters &old_parameters) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;

  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  void simulate_into(RandomNumberGenerator &rng,
                     const Particle &particle,
                     ProposalTarget &target) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;

  arma::colvec move_one_component(RandomNumberGenerator &rng,
                                  const Particle &particle) const;

  void make_copy(const ComponentwiseRandomWalkProposalKernel &another);

  std::string variable;
  arma::colvec component_scales;

};

}

#endif
