#ifndef GAUSSIANINDEPENDENTPROPOSALKERNEL_H
#define GAUSSIANINDEPENDENTPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "independent_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file gaussian_independent_proposal_kernel.h
   * @brief Defines the GaussianIndependentProposalKernel class.
   *
   * An independent gaussian proposal kernel. Proposes new parameter values independently of the current state by sampling from a gaussian distribution.
   *
   * @namespace ilike
   * @class GaussianIndependentProposalKernel
   * @brief A gaussian independent proposal kernel derived from IndependentProposalKernel.
   */


class GaussianIndependentProposalKernel : public IndependentProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for GaussianIndependentProposalKernel.
   */
  GaussianIndependentProposalKernel();
  /**
   * @brief Destructor for GaussianIndependentProposalKernel.
   */
  virtual ~GaussianIndependentProposalKernel();
  
  // find mean and cov adaptively
  /**
   * @brief Constructs a GaussianIndependentProposalKernel object.
   *
   * @param variable_names_in The variable names.
   */
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // find cov adaptively
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::colvec> &means_in);
  
  // find mean adaptively
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::mat> &covariances_in);
  
  GaussianIndependentProposalKernel(const std::string &variable_name_in,
                                    const double &mean_in,
                                    const double &sd_in);
  
  GaussianIndependentProposalKernel(const std::string &variable_name_in,
                                    const arma::colvec &mean_in,
                                    const arma::mat &covariance_in);
  
  GaussianIndependentProposalKernel(const std::vector<std::string> &variable_names_in,
                                    const std::vector<arma::colvec> &means_in,
                                    const std::vector<arma::mat> &covariances_in);
  
  /**
   * @brief Copy constructor for GaussianIndependentProposalKernel.
   *
   * @param another The GaussianIndependentProposalKernel instance to copy from.
   */
  GaussianIndependentProposalKernel(const GaussianIndependentProposalKernel &another);
  
  /**
   * @brief Assignment operator for GaussianIndependentProposalKernel.
   *
   * @param another The GaussianIndependentProposalKernel instance to copy from.
   */
  void operator=(const GaussianIndependentProposalKernel &another);
  /**
   * @brief Creates a deep copy of this GaussianIndependentProposalKernel object.
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
   * @brief Creates a deep copy and returns it as a independent_proposal_kernel pointer.
   *
   * @return The result.
   */
  IndependentProposalKernel* independent_proposal_kernel_duplicate() const;
  
  void set_mean(const std::string &variable,
                const arma::colvec &mean_in);
  
  void set_covariance(const std::string &variable,
                      const arma::mat &covariance_in);
  
  /**
   * @brief Returns the inverse covariance.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat get_inverse_covariance(const std::string &variable) const;
  /**
   * @brief Returns the covariance.
   *
   * @param variable The variable.
   *
   * @return The result.
   */
  arma::mat get_covariance(const std::string &variable) const;
  /**
   * @brief Returns the covariance.
   *
   * @param variables The variables.
   *
   * @return The result.
   */
  arma::mat get_covariance(const std::vector<std::string> &variables) const;
  
  /**
   * @brief Evaluates the independent kernel.
   *
   * @param proposed_particle The proposed particle.
   *
   * @return The result.
   */
  double evaluate_independent_kernel(const Parameters &proposed_particle) const;
  //double evaluate_independent_kernel(Variables* proposed_particle,
  //                                   const Parameters &conditioned_on_parameters) const;
  /**
   * @brief Performs the subsample evaluate independent kernel operation.
   *
   * @param proposed_particle The proposed particle.
   *
   * @return The result.
   */
  double subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const;
  
  /**
   * @brief Performs the independent simulate operation.
   *
   * @param rng The rng.
   *
   * @return The result.
   */
  Parameters independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters independent_simulate(RandomNumberGenerator &rng,
                                  const Parameters &conditioned_on_parameters) const;
  
  /**
   * @brief Performs the subsample independent simulate operation.
   *
   * @param rng The rng.
   *
   * @return The result.
   */
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable) const;
  
  Parameters subsample_independent_simulate(RandomNumberGenerator &rng,
                                            const std::string &variable,
                                            const Parameters &conditioned_on_parameters) const;
  
  arma::mat independent_gradient_of_log(const std::string &variable,
                                        const Parameters &proposed_particle);
  arma::mat subsample_independent_gradient_of_log(const std::string &variable,
                                                  const Parameters &proposed_particle);
  
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
  // MH has sim prop and eval prop, t**ake in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  /**
   * @brief Copies the state of another GaussianIndependentProposalKernel into this object.
   *
   * @param another The GaussianIndependentProposalKernel instance to copy from.
   */
  void make_copy(const GaussianIndependentProposalKernel &another);
  
  boost::unordered_map< std::string, GaussianProposalInfo> proposal_info;
  
};
}

#endif
