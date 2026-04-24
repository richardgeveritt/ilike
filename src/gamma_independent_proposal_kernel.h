#ifndef GAMMAINDEPENDENTPROPOSALKERNEL_H
#define GAMMAINDEPENDENTPROPOSALKERNEL_H

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
   * @file gamma_independent_proposal_kernel.h
   * @brief Defines the GammaIndependentProposalKernel class.
   *
   * An independent gamma proposal kernel. Proposes new parameter values independently of the current state by sampling from a gamma distribution.
   *
   * @namespace ilike
   * @class GammaIndependentProposalKernel
   * @brief A gamma independent proposal kernel derived from IndependentProposalKernel.
   */


class GammaIndependentProposalKernel : public IndependentProposalKernel
{
  
public:
  
  /**
   * @brief Default constructor for GammaIndependentProposalKernel.
   */
  GammaIndependentProposalKernel();
  /**
   * @brief Destructor for GammaIndependentProposalKernel.
   */
  virtual ~GammaIndependentProposalKernel();
  
  GammaIndependentProposalKernel(const std::string &variable_in,
                                 const double &shape_in,
                                 const double &rate_in);
  
  /**
   * @brief Copy constructor for GammaIndependentProposalKernel.
   *
   * @param another The GammaIndependentProposalKernel instance to copy from.
   */
  GammaIndependentProposalKernel(const GammaIndependentProposalKernel &another);
  
  /**
   * @brief Assignment operator for GammaIndependentProposalKernel.
   *
   * @param another The GammaIndependentProposalKernel instance to copy from.
   */
  void operator=(const GammaIndependentProposalKernel &another);
  /**
   * @brief Creates a deep copy of this GammaIndependentProposalKernel object.
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
   * @brief Copies the state of another GammaIndependentProposalKernel into this object.
   *
   * @param another The GammaIndependentProposalKernel instance to copy from.
   */
  void make_copy(const GammaIndependentProposalKernel &another);
  
  /** @brief The variable. */
  std::string variable;
  /** @brief The shape. */
  double shape;
  /** @brief The rate. */
  double rate;
  
};
}

#endif
