#ifndef PRINCIPALCOMPONENTSPROPOSALKERNEL_H
#define PRINCIPALCOMPONENTSPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "gaussian_random_walk_proposal_kernel.h"
#include "principal_components_info.h"

namespace ilike
{

// Principal-components Metropolis update (Andrieu & Thoms 2008, Algorithm 8;
// dev/plan.md 3.2, A4): a direction l, chosen uniformly among the first m
// principal directions of Sigma, and a move N(0, ell_l rho_l) along it. The
// directions start as the eigenvectors of the covariance given, with all
// ell_l equal to the scale given. Symmetric; its density is singular (used in
// m_proposal only, where it cancels).
class PrincipalComponentsProposalKernel : public GaussianRandomWalkProposalKernel
{

public:

  PrincipalComponentsProposalKernel();
  virtual ~PrincipalComponentsProposalKernel();

  // number_of_directions 0: all d directions
  PrincipalComponentsProposalKernel(const std::string &variable_name_in,
                                    const arma::mat &covariance_in,
                                    double scale_in,
                                    size_t number_of_directions);

  PrincipalComponentsProposalKernel(const PrincipalComponentsProposalKernel &another);
  void operator=(const PrincipalComponentsProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;

  PrincipalComponentsInfo* get_principal_components(const std::string &variable);

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

  arma::colvec move_along_one_direction(RandomNumberGenerator &rng,
                                        const Particle &particle) const;

  void make_copy(const PrincipalComponentsProposalKernel &another);

  std::string variable;
  PrincipalComponentsInfo components;

};

}

#endif
