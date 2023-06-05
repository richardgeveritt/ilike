#ifndef COMPOSITEPROPOSALKERNEL_H
#define COMPOSITEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"

class CompositeProposalKernel : public ProposalKernel
{

public:

  CompositeProposalKernel();
  virtual ~CompositeProposalKernel();

  CompositeProposalKernel(const std::vector<ProposalKernel*> &all_kernels_in);
  
  CompositeProposalKernel(const CompositeProposalKernel &another);

  void operator=(const CompositeProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:
  
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle) const;
  
  /*
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle,
                                  const Parameters &conditioned_on_parameters) const;
  */
  
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle) const;
  
  /*
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle,
                                            const Parameters &conditioned_on_parameters) const;
  */
  
  Parameters simulate(RandomNumberGenerator &rng,
                      Particle &particle) const;
  
  /*
  Parameters simulate(RandomNumberGenerator &rng,
                      Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                Particle &particle) const;
  
  /*
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                Particle &particle) const;
  
  /*
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  */
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle);
  
  /*
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle,
                                     const Parameters &conditioned_on_parameters);
  */

  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;
  
  /*
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle,
                                               const Parameters &conditioned_on_parameters);
  */
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle);

  void make_copy(const CompositeProposalKernel &another);
  
  // stored here
  std::vector<ProposalKernel*> all_kernels;
  
};

#endif
