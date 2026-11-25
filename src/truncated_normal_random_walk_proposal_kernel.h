#ifndef TRUNCATEDNORMALRANDOMWALKPROPOSALKERNEL_H
#define TRUNCATEDNORMALRANDOMWALKPROPOSALKERNEL_H

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

// The random walk tnorm_rw(x, sd, lower, upper) (C217): each element moves
// to a draw from the normal(x, sd) truncated to [lower, upper], so a
// proposal never leaves the limits. It is not symmetric: the truncated mass
// Phi((upper - x)/sd) - Phi((lower - x)/sd) depends on the current point, so
// q(x* | x) is in the Metropolis-Hastings ratio (an MH proposal only).
class TruncatedNormalRandomWalkProposalKernel : public ProposalKernel
{
public:

  TruncatedNormalRandomWalkProposalKernel();
  virtual ~TruncatedNormalRandomWalkProposalKernel();
  TruncatedNormalRandomWalkProposalKernel(const std::string &variable_in, double sd_in,
                                          double lower_in, double upper_in);
  TruncatedNormalRandomWalkProposalKernel(const TruncatedNormalRandomWalkProposalKernel &another);
  void operator=(const TruncatedNormalRandomWalkProposalKernel &another);
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

  void make_copy(const TruncatedNormalRandomWalkProposalKernel &another);
  std::string variable;
  double sd;
  double lower;
  double upper;
};

}

#endif
