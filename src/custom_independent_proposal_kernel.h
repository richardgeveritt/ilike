#ifndef CUSTOMINDEPENDENTPROPOSALKERNEL_H
#define CUSTOMINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "function_pointers.h"

class CustomIndependentProposalKernel : public IndependentProposalKernel
{

public:

  CustomIndependentProposalKernel();
  virtual ~CustomIndependentProposalKernel();
  
  CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in);
  
  CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
                                  EvaluateLogIndependentProposalPtr proposal_evaluate_in);
  
  CustomIndependentProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
                                  EvaluateLogIndependentProposalPtr proposal_evaluate_in,
                                  const Parameters &proposal_parameters_in);

  CustomIndependentProposalKernel(const CustomIndependentProposalKernel &another);

  void operator=(const CustomIndependentProposalKernel &another);
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
  //arma::mat independent_gradient_of_log(const std::string &variable,
  //                                      Variables* proposed_particle,
  //                                      const Parameters &conditioned_on_parameters);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  //arma::mat subsample_independent_gradient_of_log(const std::string &variable,
  //                                                Variables* proposed_particle,
  //                                                const Parameters &conditioned_on_parameters);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:

  void make_copy(const CustomIndependentProposalKernel &another);
  
  EvaluateLogIndependentProposalPtr proposal_evaluate;
  
  SimulateIndependentProposalPtr proposal_simulate;
  
  Parameters proposal_parameters;
  
};

#endif