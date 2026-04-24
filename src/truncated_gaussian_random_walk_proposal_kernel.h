#ifndef TRUNCATEDGAUSSIANRANDOMWALKPROPOSALKERNEL_H
#define TRUNCATEDGAUSSIANRANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <boost/unordered_map.hpp>
#include <vector>

#include "distributions.h"
#include "gaussian_proposal_info.h"
#include "ilike_header.h"
#include "particle.h"
#include "symmetric_proposal_kernel.h"

namespace ilike {
  /**
   * @file truncated_gaussian_random_walk_proposal_kernel.h
   * @brief Defines the TruncatedGaussianRandomWalkProposalKernel class.
   *
   * A truncated gaussian random walk proposal kernel. Proposes new parameter values during MCMC or SMC moves using a truncated gaussian random walk distribution centred on the current state.
   *
   * @namespace ilike
   * @class TruncatedGaussianRandomWalkProposalKernel
   * @brief A truncated gaussian random walk proposal kernel derived from SymmetricProposalKernel.
   */


class TruncatedGaussianRandomWalkProposalKernel
    : public SymmetricProposalKernel {

public:
  /**
   * @brief Default constructor for TruncatedGaussianRandomWalkProposalKernel.
   */
  TruncatedGaussianRandomWalkProposalKernel();
  /**
   * @brief Destructor for TruncatedGaussianRandomWalkProposalKernel.
   */
  virtual ~TruncatedGaussianRandomWalkProposalKernel();

  // make cov_names from var_names and find cov adaptively
  TruncatedGaussianRandomWalkProposalKernel(
      const std::vector<std::string> &variable_names_in,
      const std::vector<arma::colvec> &lower_in,
      const std::vector<arma::colvec> &upper_in);

  // make cov_names from var_names
  TruncatedGaussianRandomWalkProposalKernel(
      const std::vector<std::string> &variable_names_in,
      const std::vector<arma::mat> &covariances_in,
      const std::vector<arma::colvec> &lower_in,
      const std::vector<arma::colvec> &upper_in);

  // make cov_names from var_names
  TruncatedGaussianRandomWalkProposalKernel(const std::string &variable_name_in,
                                            const arma::mat &covariance_in,
                                            const arma::colvec &lower_in,
                                            const arma::colvec &upper_in);

  TruncatedGaussianRandomWalkProposalKernel(const std::string &variable_name_in,
                                            const arma::mat &covariance_in,
                                            double scale_in,
                                            const arma::colvec &lower_in,
                                            const arma::colvec &upper_in);

  TruncatedGaussianRandomWalkProposalKernel(const std::string &variable_name_in,
                                            const double &sd_in,
                                            const double &lower_in,
                                            const double &upper_in);

  TruncatedGaussianRandomWalkProposalKernel(
      const TruncatedGaussianRandomWalkProposalKernel &another);

  /**
   * @brief Assignment operator for TruncatedGaussianRandomWalkProposalKernel.
   *
   * @param another The TruncatedGaussianRandomWalkProposalKernel instance to copy from.
   */
  void operator=(const TruncatedGaussianRandomWalkProposalKernel &another);
  Kernel *duplicate() const;
  ProposalKernel *proposal_kernel_duplicate() const;
  SymmetricProposalKernel *symmetric_proposal_kernel_duplicate() const;

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
  void set_proposal_parameters(Parameters *proposal_parameters_in);

  /**
   * @brief Returns the variables.
   *
   * @return The result.
   */
  std::vector<std::string> get_variables() const;

  GradientEstimatorOutput *simulate_gradient_estimator_output() const;

  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  std::vector<const ProposalKernel *> get_proposals() const;

  /**
   * @brief Sets the index.
   *
   * @param index_in The index.
   */
  void set_index(Index *index_in);
  /**
   * @brief Sets the index if null.
   *
   * @param index_in The index.
   */
  void set_index_if_null(Index *index_in);

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
  void set_data(Data *data_in);

  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in
  // acceptance, Set current value if accepted. Proposal needs to call simulate
  // in all llhdoutputs

protected:
  // bool unused_variables_kept;

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

  arma::mat
  specific_subsample_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);

  /*
   arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */

  void make_copy(const TruncatedGaussianRandomWalkProposalKernel &another);

  boost::unordered_map<std::string, GaussianProposalInfo>
      proposal_info;

  boost::unordered_map<std::string, arma::colvec> lower_info;
  boost::unordered_map<std::string, arma::colvec> upper_info;
};
} // namespace ilike

#endif
