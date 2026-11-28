#ifndef DIMENSIONJUMPPROPOSALKERNEL_H
#define DIMENSIONJUMPPROPOSALKERNEL_H

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

// (dev/plan.md 3.22 U6, C313) The proposal of a reversible-jump move on a
// model indicator k in {minimum, ..., maximum}: k + 1 or k - 1 with
// probability 1/2 each, k + 1 at the minimum and k - 1 at the maximum. Its
// density ratio q(k | k')/q(k' | k) (1/2 against 1 at the bounds) is in the
// Metropolis-Hastings ratio; the move's step makes the dimension matching.
class DimensionJumpProposalKernel : public ProposalKernel
{
public:

  DimensionJumpProposalKernel();
  virtual ~DimensionJumpProposalKernel();
  DimensionJumpProposalKernel(const std::string &variable_in, double minimum_in, double maximum_in);
  DimensionJumpProposalKernel(const DimensionJumpProposalKernel &another);
  void operator=(const DimensionJumpProposalKernel &another);
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

  void make_copy(const DimensionJumpProposalKernel &another);
  std::string variable;
  double minimum;
  double maximum;
  // log q(to | from)
  double log_density(double to, double from) const;
};

}

#endif
