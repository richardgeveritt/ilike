#ifndef LANGEVINPROPOSALKERNEL_H
#define LANGEVINPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class GradientEstimator;
class Index;

class LangevinProposalKernel : public ProposalKernel
{

public:

  LangevinProposalKernel();
  virtual ~LangevinProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  LangevinProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // make cov_names from var_names
  LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                         const std::vector<arma::mat> &covariances_in);
  
  LangevinProposalKernel(const std::string &variable_name_in,
                         const arma::mat &covariance_in);

  LangevinProposalKernel(const LangevinProposalKernel &another);

  void operator=(const LangevinProposalKernel &another);
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
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle);
  
  /*
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle,
                                               const Parameters &conditioned_on_parameters);
  */

  void make_copy(const LangevinProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
  bool unused_variables_kept;
  
  // stored here
  GradientEstimator* gradient_estimator;
  
  // stored here
  Index* index;

};

#endif
