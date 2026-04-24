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

namespace ilike
{
  /**
   * @file langevin_proposal_kernel.h
   * @brief Defines the GradientEstimator class.
   *
   * Estimates the gradient for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class GradientEstimator
   * @brief The gradient estimator class.
   */



class GradientEstimator;
class Index;

class LangevinProposalKernel : public ProposalKernel
{
  
public:
  
  /**
   * @brief Performs the langevinproposalkernel operation.
   */
  LangevinProposalKernel();
  /**
   * @brief Performs the ~langevinproposalkernel operation.
   */
  virtual ~LangevinProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                         GradientEstimator* gradient_estimator_in);
  
  // make cov_names from var_names
  LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                         const std::vector<arma::mat> &covariances_in,
                         GradientEstimator* gradient_estimator_in);
  
  LangevinProposalKernel(const std::string &variable_name_in,
                         const arma::mat &covariance_in,
                         GradientEstimator* gradient_estimator_in);
  
  LangevinProposalKernel(const std::string &variable_name_in,
                         const arma::mat &covariance_in,
                         double scale_in,
                         GradientEstimator* gradient_estimator_in);
  
  LangevinProposalKernel(const std::string &variable_name_in,
                         const double &sd_in,
                         GradientEstimator* gradient_estimator_in);
  
  /**
   * @brief Performs the langevinproposalkernel operation.
   *
   * @param another The GradientEstimator instance to copy from.
   */
  LangevinProposalKernel(const LangevinProposalKernel &another);
  
  /**
   * @brief Assignment operator for GradientEstimator.
   *
   * @param another The GradientEstimator instance to copy from.
   */
  void operator=(const LangevinProposalKernel &another);
  /**
   * @brief Creates a deep copy of this GradientEstimator object.
   *
   * @return The result.
   */
  Kernel* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a proposal_kernel pointer.
   *
   * @return The result.
   */
  ProposalKernel* proposal_kernel_duplicate() const;
  
  /**
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  /**
   * @brief Simulates gradient estimator output.
   *
   * @return The result.
   */
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  
  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  std::vector<const ProposalKernel*> get_proposals() const;
  
  /**
   * @brief Sets the index.
   *
   * @param index_in The index.
   */
  void set_index(Index* index_in);
  /**
   * @brief Sets the index if null.
   *
   * @param index_in The index.
   */
  void set_index_if_null(Index* index_in);
  
  /**
   * @brief Performs the can be evaluated operation.
   *
   * @return The result.
   */
  bool can_be_evaluated() const;
  
  /**
   * @brief Sets the data.
   *
   * @param data_in The data.
   */
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
  
  void make_copy(const LangevinProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
  //bool unused_variables_kept;
  
  // stored here
  /** @brief The gradient estimator. */
  GradientEstimator* gradient_estimator;
  
  // stored here
  /** @brief The index. */
  Index* index;
  
};
}

#endif
