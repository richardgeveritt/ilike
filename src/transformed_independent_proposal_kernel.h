#ifndef TRANSFORMEDINDEPENDENTPROPOSALKERNEL_H
#define TRANSFORMEDINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "transform.h"

class TransformedIndependentProposalKernel : public IndependentProposalKernel
{

public:

  TransformedIndependentProposalKernel();
  virtual ~TransformedIndependentProposalKernel();
  
  TransformedIndependentProposalKernel(IndependentProposalKernel* proposal_in,
                                       const Transform &transform_in,
                                       bool distribution_on_transformed_space_in=true);

  TransformedIndependentProposalKernel(const TransformedIndependentProposalKernel &another);

  void operator=(const TransformedIndependentProposalKernel &another);
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
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, t**ake in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const TransformedIndependentProposalKernel &another);
  
  // stored here
  IndependentProposalKernel* proposal;
  
  Transform transform;
  
  bool distribution_on_transformed_space;
  
};

#endif
