#ifndef MIRRORPROPOSALKERNEL_H
#define MIRRORPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

class MirrorProposalKernel : public SymmetricProposalKernel
{

public:

  MirrorProposalKernel();
  virtual ~MirrorProposalKernel();
  
  MirrorProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // find cov adaptively
  MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                       const std::vector<arma::colvec> &means_in);
  
  // find mean adaptively
  MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                       const std::vector<arma::mat> &covariances_in);
  
  MirrorProposalKernel(const std::string &variable_name_in,
                       const arma::colvec &mean_in,
                       const arma::mat &covariance_in);
  
  MirrorProposalKernel(const std::vector<std::string> &variable_names_in,
                       const std::vector<arma::colvec> &means_in,
                       const std::vector<arma::mat> &covariances_in);

  MirrorProposalKernel(const MirrorProposalKernel &another);

  void operator=(const MirrorProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
  arma::mat get_inverse_covariance(const std::string &variable);
  
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

  void make_copy(const MirrorProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;

};

#endif
