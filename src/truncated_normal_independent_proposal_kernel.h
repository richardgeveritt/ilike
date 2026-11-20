#ifndef TRUNCATEDNORMALINDEPENDENTPROPOSALKERNEL_H
#define TRUNCATEDNORMALINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

// The independent proposal tnorm(x, mean, sd, lower, upper) (C217): draws
// from the normal(mean, sd) truncated to [lower, upper] (Robert's sampler,
// distributions.h rtnorm), for one value or each element of a block.
class TruncatedNormalIndependentProposalKernel : public IndependentProposalKernel
{
public:

  TruncatedNormalIndependentProposalKernel();
  virtual ~TruncatedNormalIndependentProposalKernel();
  TruncatedNormalIndependentProposalKernel(const std::string &variable_in, double mean_in, double sd_in,
                                           double lower_in, double upper_in);
  TruncatedNormalIndependentProposalKernel(const TruncatedNormalIndependentProposalKernel &another);
  void operator=(const TruncatedNormalIndependentProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;

  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  Parameters independent_simulate_like(RandomNumberGenerator &rng,
                                       const Parameters &shapes) const;
  Parameters independent_simulate(RandomNumberGenerator &rng,
                                  const Parameters &conditioned_on_parameters) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const Parameters &conditioned_on_parameters) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable,
                                            const Parameters &conditioned_on_parameters) const;
  arma::mat independent_gradient_of_log(const std::string &variable,
                                        const Parameters &proposed_particle);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

protected:

  void make_copy(const TruncatedNormalIndependentProposalKernel &another);
  arma::colvec draw(RandomNumberGenerator &rng, arma::uword n) const;
  std::string variable;
  double mean;
  double sd;
  double lower;
  double upper;
};

}

#endif
