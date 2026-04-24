#ifndef NONLINEARGAUSSIANNOISEPROPOSALKERNEL_H
#define NONLINEARGAUSSIANNOISEPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "gaussian_noise_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"
#include "transform.h"

namespace ilike
{
  /**
   * @file nonlinear_gaussian_noise_proposal_kernel.h
   * @brief Defines the NonLinearGaussianNoiseProposalKernel class.
   *
   * A non linear gaussian noise proposal kernel. Proposes new parameter values during MCMC or SMC moves using a non linear gaussian noise distribution centred on the current state.
   *
   * @namespace ilike
   * @class NonLinearGaussianNoiseProposalKernel
   * @brief A non linear gaussian noise proposal kernel derived from GaussianNoiseProposalKernel.
   */


class NonLinearGaussianNoiseProposalKernel : public GaussianNoiseProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for NonLinearGaussianNoiseProposalKernel.
   */
  NonLinearGaussianNoiseProposalKernel();
  /**
   * @brief Destructor for NonLinearGaussianNoiseProposalKernel.
   */
  virtual ~NonLinearGaussianNoiseProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  NonLinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                       std::shared_ptr<Transform> transform_in);
  
  // make cov_names from var_names
  NonLinearGaussianNoiseProposalKernel(const std::vector<std::string> &variable_names_in,
                                       std::shared_ptr<Transform> transform_in,
                                       const std::vector<arma::mat> &covariances_in);
  
  // make cov_names from var_names
  NonLinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                       std::shared_ptr<Transform> transform_in,
                                       const arma::mat &covariance_in);
  
  NonLinearGaussianNoiseProposalKernel(const std::string &variable_name_in,
                                       std::shared_ptr<Transform> transform_in,
                                       const double &sd_in);
  
  /**
   * @brief Copy constructor for NonLinearGaussianNoiseProposalKernel.
   *
   * @param another The NonLinearGaussianNoiseProposalKernel instance to copy from.
   */
  NonLinearGaussianNoiseProposalKernel(const NonLinearGaussianNoiseProposalKernel &another);
  
  /**
   * @brief Assignment operator for NonLinearGaussianNoiseProposalKernel.
   *
   * @param another The NonLinearGaussianNoiseProposalKernel instance to copy from.
   */
  void operator=(const NonLinearGaussianNoiseProposalKernel &another);
  /**
   * @brief Creates a deep copy of this NonLinearGaussianNoiseProposalKernel object.
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
   * @brief Creates a deep copy and returns it as a gaussian_noise_proposal_kernel pointer.
   *
   * @return The result.
   */
  GaussianNoiseProposalKernel* gaussian_noise_proposal_kernel_duplicate() const;
  
  void set_covariance(const std::string &variable,
                      const arma::mat &covariance_in);
  
  /**
   * @brief Returns the inverse covariance.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat get_inverse_covariance(const std::string &variable);
  /**
   * @brief Returns the covariance.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat get_covariance(const std::string &variable);
  
  /**
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  void set_proposal_parameters(Parameters* proposal_parameters_in);
  
  /**
   * @brief Returns the variables.
   *
   * @return The result.
   */
  std::vector<std::string> get_variables() const;
  
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
  
  void make_copy(const NonLinearGaussianNoiseProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  /** @brief The transform. */
  std::shared_ptr<Transform> transform;
  
};
}

#endif
