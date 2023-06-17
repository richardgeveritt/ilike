#ifndef CUSTOMSYMMETRICPROPOSALKERNEL_H
#define CUSTOMSYMMETRICPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

class CustomSymmetricProposalKernel : public SymmetricProposalKernel
{

public:

  CustomSymmetricProposalKernel();
  virtual ~CustomSymmetricProposalKernel();
  
  CustomSymmetricProposalKernel(SimulateMCMCProposalPtr proposal_simulate_in);

  CustomSymmetricProposalKernel(const CustomSymmetricProposalKernel &another);

  void operator=(const CustomSymmetricProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:
  
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle) const;
  
  
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle) const;
  
  Parameters simulate(RandomNumberGenerator &rng,
                      Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                Particle &particle) const;
  
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle);
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle);
  
  
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;

  void make_copy(const CustomSymmetricProposalKernel &another);
  
  EvaluateLogMCMCProposalPtr proposal_evaluate;
  
  SimulateMCMCProposalPtr proposal_simulate;
  
  Parameters* proposal_parameters;
  
  
};

#endif