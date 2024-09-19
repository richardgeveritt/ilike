#ifndef LINEARGAUSSIANNOISEPROPOSALKERNEL_H
#define LINEARGAUSSIANNOISEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "gaussian_noise_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
class LinearGaussianNoiseProposalKernel : public GaussianNoiseProposalKernel
{
  
public:
  
  LinearGaussianNoiseProposalKernel();
  virtual ~LinearGaussianNoiseProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  LinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<std::string> &conditioned_on_variable_names_in,
                                    const std::vector<arma::mat> &As_in);
  
  // make cov_names from var_names
  LinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<std::string> &conditioned_on_variable_names_in,
                                    const std::vector<arma::mat> &As_in,
                                    const std::vector<arma::mat> &covariances_in);
  
  // make cov_names from var_names
  LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                    const std::vector<std::string> &conditioned_on_variable_names_in,
                                    const arma::mat &A_in,
                                    const arma::mat &covariance_in);
  
  LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                    const std::string &conditioned_on_variable_name_in,
                                    const arma::mat &A_in,
                                    const arma::mat &covariance_in);
  
  LinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                    const std::vector<std::string> &conditioned_on_variable_names_in,
                                    const double &sd_in);
  
  LinearGaussianNoiseProposalKernel(const LinearGaussianNoiseProposalKernel &another);
  
  void operator=(const LinearGaussianNoiseProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;
  GaussianNoiseProposalKernel* gaussian_noise_proposal_kernel_duplicate() const;
  
  void set_covariance(const std::string &variable,
                      const arma::mat &covariance_in);
  
  arma::mat get_inverse_covariance(const std::string &variable);
  arma::mat get_covariance(const std::string &variable);
  
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  std::vector<std::string> get_variables() const;
  
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  
  std::vector<const ProposalKernel*> get_proposals() const;
  
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  //bool unused_variables_kept;
  
  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  
  /*
   double specific_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  
  /*
   double specific_subsample_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  
  /*
   Parameters simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  
  /*
   Parameters subsample_simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  
  /*
   Parameters subsample_simulate(RandomNumberGenerator &rng,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  /*
   arma::mat specific_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  /*
   arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void make_copy(const LinearGaussianNoiseProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  std::vector<std::string> conditioned_on_variable_names;
};
}

#endif
