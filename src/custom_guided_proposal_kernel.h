#ifndef CUSTOMGUIDEDPROPOSALKERNEL_H
#define CUSTOMGUIDEDPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
class CustomGuidedProposalKernel : public ProposalKernel
{

public:

  CustomGuidedProposalKernel();
  virtual ~CustomGuidedProposalKernel();
  
  CustomGuidedProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                             Data* data_in);
  
  CustomGuidedProposalKernel(SimulateGuidedMCMCProposalPtr proposal_simulate_in,
                             EvaluateLogGuidedMCMCProposalPtr proposal_evaluate_in,
                             Data* data_in);

  CustomGuidedProposalKernel(const CustomGuidedProposalKernel &another);

  void operator=(const CustomGuidedProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  
  bool can_be_evaluated() const;
  
  void set_data(Data* data_in);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

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
  
  
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;

  void make_copy(const CustomGuidedProposalKernel &another);
  
  EvaluateLogGuidedMCMCProposalPtr proposal_evaluate;
  
  SimulateGuidedMCMCProposalPtr proposal_simulate;
  
  Parameters* proposal_parameters;
  
  Data* data;
};
}

#endif
