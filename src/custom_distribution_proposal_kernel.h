#ifndef CUSTOMDISTRIBUTIONPROPOSALKERNEL_H
#define CUSTOMDISTRIBUTIONPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "function_pointers.h"

class CustomDistributionProposalKernel : public IndependentProposalKernel
{

public:

  CustomDistributionProposalKernel();
  virtual ~CustomDistributionProposalKernel();
  
  CustomDistributionProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in);
  
  CustomDistributionProposalKernel(SimulateIndependentProposalPtr proposal_simulate_in,
                                   EvaluateLogDistributionPtr proposal_evaluate_in);

  CustomDistributionProposalKernel(const CustomDistributionProposalKernel &another);

  void operator=(const CustomDistributionProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;
  
  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double evaluate_independent_kernel(Variables* proposed_particle,
  //                                   const Parameters &conditioned_on_parameters) const;
  
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double subsample_evaluate_independent_kernel(Variables* proposed_particle,
  //                                             const Parameters &conditioned_on_parameters) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  
  //Parameters independent_simulate(RandomNumberGenerator &rng,
  //                                const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const;
  
  //Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
  //                                          const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng,
                                  const Parameters &conditioned_on_parameters) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const Parameters &conditioned_on_parameters) const;
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable,
                                            const Parameters &conditioned_on_parameters) const;  
  
  //Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
  //                                          const std::string &variable,
  //                                          const Parameters &conditioned_on_parameters) const;
  
  /*
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Variables* proposed_particle,
                                     Variables* old_particle);
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Variables* proposed_particle,
                                     Variables* old_particle,
                                     const Parameters &conditioned_on_parameters);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Variables* proposed_particle,
                                               Variables* old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Variables* proposed_particle,
                                               Variables* old_particle,
                                               const Parameters &conditioned_on_parameters);
  */
  
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

  void make_copy(const CustomDistributionProposalKernel &another);
  
  EvaluateLogDistributionPtr proposal_evaluate;
  
  SimulateIndependentProposalPtr proposal_simulate;
  
  Parameters proposal_parameters;
  
};

#endif