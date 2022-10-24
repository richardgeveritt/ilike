#ifndef GAUSSIANRANDOMWALKPROPOSALKERNEL_H
#define GAUSSIANRANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "function_pointers.h"
#include "gaussian_proposal_info.h"

class GaussianRandomWalkProposalKernel : public SymmetricProposalKernel
{

public:

  GaussianRandomWalkProposalKernel();
  virtual ~GaussianRandomWalkProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  GaussianRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // make cov_names from var_names
  GaussianRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                   const std::vector<arma::mat> &covariances_in);

  GaussianRandomWalkProposalKernel(const GaussianRandomWalkProposalKernel &another);

  void operator=(const GaussianRandomWalkProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
  void set_covariance(const std::string &variable,
                      const arma::mat &covariance_in);
  
  arma::mat get_inverse_covariance(const std::string &variable);
  arma::mat get_covariance(const std::string &variable);
  
// Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs

protected:
  
  bool unused_variables_kept;
  
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle) const;
  double specific_evaluate_kernel(Particle &proposed_particle,
                                  Particle &old_particle,
                                  const Parameters &conditioned_on_parameters) const;
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                            Particle &old_particle,
                                            const Parameters &conditioned_on_parameters) const;
  
  Parameters simulate(RandomNumberGenerator &rng,
                      Particle &particle) const;
  
  Parameters simulate(RandomNumberGenerator &rng,
                      Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                      Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle);
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     Particle &proposed_particle,
                                     Particle &old_particle,
                                     const Parameters &conditioned_on_parameters);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               Particle &proposed_particle,
                                               Particle &old_particle,
                                               const Parameters &conditioned_on_parameters);

  void make_copy(const GaussianRandomWalkProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;

};

#endif