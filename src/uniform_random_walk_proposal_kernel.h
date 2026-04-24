#ifndef UNIFORMRANDOMWALKPROPOSALKERNEL_H
#define UNIFORMRANDOMWALKPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <boost/unordered_map.hpp>

#include "symmetric_proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file uniform_random_walk_proposal_kernel.h
   * @brief Defines the UniformRandomWalkProposalKernel class.
   *
   * A uniform random walk proposal kernel. Proposes new parameter values during MCMC or SMC moves using a uniform random walk distribution centred on the current state.
   *
   * @namespace ilike
   * @class UniformRandomWalkProposalKernel
   * @brief An uniform random walk proposal kernel derived from SymmetricProposalKernel.
   */


class UniformRandomWalkProposalKernel : public SymmetricProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for UniformRandomWalkProposalKernel.
   */
  UniformRandomWalkProposalKernel();
  /**
   * @brief Destructor for UniformRandomWalkProposalKernel.
   */
  virtual ~UniformRandomWalkProposalKernel();
  
  // make cov_names from var_names and find cov adaptively
  /**
   * @brief Constructs a UniformRandomWalkProposalKernel object.
   *
   * @param variable_names_in The variable names.
   */
  UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in);
  
  // make cov_names from var_names
  UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                  const std::vector<arma::mat> &halfwidths_in);
  
  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const double &halfwidth_in);
  
  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const arma::colvec &halfwidth_in);
  
  UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                  const arma::mat &halfwidth_in);
  
  /**
   * @brief Copy constructor for UniformRandomWalkProposalKernel.
   *
   * @param another The UniformRandomWalkProposalKernel instance to copy from.
   */
  UniformRandomWalkProposalKernel(const UniformRandomWalkProposalKernel &another);
  
  /**
   * @brief Assignment operator for UniformRandomWalkProposalKernel.
   *
   * @param another The UniformRandomWalkProposalKernel instance to copy from.
   */
  void operator=(const UniformRandomWalkProposalKernel &another);
  /**
   * @brief Creates a deep copy of this UniformRandomWalkProposalKernel object.
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
   * @brief Creates a deep copy and returns it as a symmetric_proposal_kernel pointer.
   *
   * @return The result.
   */
  SymmetricProposalKernel* symmetric_proposal_kernel_duplicate() const;
  
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
  
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle,
                      const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle,
                                const Parameters &conditioned_on_parameters) const;
  
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
  
  void make_copy(const UniformRandomWalkProposalKernel &another);
  
  //bool unused_variables_kept;
  
  // store the half widths for every variable
  boost::unordered_map< std::string, arma::mat> proposal_info;
  
};
}

#endif
