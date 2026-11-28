#ifndef INVOLUTIONPROPOSALKERNEL_H
#define INVOLUTIONPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "proposal_kernel.h"

namespace ilike
{

// The proposal of an involution(…) move (dev/plan.md 3.22, C316): it changes
// nothing (its density is 1 both ways). The move is a ratio-estimating MH
// move whose state step, InvolutionStep, draws the auxiliary variables and
// applies the map; the proposed particle this proposal makes is the base of
// the new state (its factors have new outputs).
class InvolutionProposalKernel : public ProposalKernel
{
public:

  InvolutionProposalKernel();
  virtual ~InvolutionProposalKernel();
  InvolutionProposalKernel(const InvolutionProposalKernel &another);
  void operator=(const InvolutionProposalKernel &another);
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

};

}

#endif
