#ifndef COMPOSITEINDEPENDENTPROPOSALKERNEL_H
#define COMPOSITEINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class CompositeIndependentProposalKernel : public IndependentProposalKernel
{

public:

  CompositeIndependentProposalKernel();
  virtual ~CompositeIndependentProposalKernel();
  
  CompositeIndependentProposalKernel(const std::vector<IndependentProposalKernel*> &all_kernels_in);

  CompositeIndependentProposalKernel(const CompositeIndependentProposalKernel &another);

  void operator=(const CompositeIndependentProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;
  
  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double evaluate_independent_kernel(Variables* proposed_particle,
  //                                   const Parameters &conditioned_on_parameters) const;
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  
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
  
  std::vector<ProposalKernel*> get_proposals();
  
  void set_index(Index* index_in);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, t**ake in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const CompositeIndependentProposalKernel &another);
  
  // stored here
  std::vector<IndependentProposalKernel*> all_kernels;
  
};

#endif
