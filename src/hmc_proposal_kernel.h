#ifndef HMCPROPOSALKERNEL_H
#define HMCPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
class HMCProposalKernel : public ProposalKernel
{
  
public:
  
  HMCProposalKernel();
  virtual ~HMCProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                    GradientEstimator* gradient_estimator_in);
  
  // make cov_names from var_names
  HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                    const std::vector<arma::mat> &covariances_in,
                    GradientEstimator* gradient_estimator_in);
  
  // find covariance adaptively
  HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                    const std::vector<std::string> &covariance_names_in,
                    GradientEstimator* gradient_estimator_in);
  
  HMCProposalKernel(const std::string &variable_name_in,
                    const arma::mat &covariance_in,
                    GradientEstimator* gradient_estimator_in);
  
  HMCProposalKernel(const std::string &variable_name_in,
                    const double &sd_in,
                    GradientEstimator* gradient_estimator_in);
  
  HMCProposalKernel(const std::vector<std::string> &variable_names_in,
                    const std::vector<std::string> &covariance_names_in,
                    const std::vector<arma::mat> &covariances_in,
                    GradientEstimator* gradient_estimator_in);
  
  HMCProposalKernel(const HMCProposalKernel &another);
  
  void operator=(const HMCProposalKernel &another);
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
  //virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
  //                                                     Particle &proposed_particle,
  //                                                     Particle &old_particle)=0;
  
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);
  
  /*
   arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  void make_copy(const HMCProposalKernel &another);
  
  //bool unused_variables_kept;
  
  std::vector<std::string> variable_names;
  std::vector<std::string> covariance_names;
  
};
}

#endif
