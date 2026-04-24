#ifndef BARKERDYNAMICSPROPOSALKERNEL_H
#define BARKERDYNAMICSPROPOSALKERNEL_H

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
   * @file barker_dynamics_proposal_kernel.h
   * @brief Defines the GradientEstimator class.
   *
   * Estimates the gradient for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class GradientEstimator
   * @brief The gradient estimator class.
   */


  class GradientEstimator;
  class IndependentProposalKernel;

  class BarkerDynamicsProposalKernel : public ProposalKernel
  {

  public:
    /**
     * @brief Performs the barkerdynamicsproposalkernel operation.
     */
    BarkerDynamicsProposalKernel();
    /**
     * @brief Performs the ~barkerdynamicsproposalkernel operation.
     */
    virtual ~BarkerDynamicsProposalKernel();

    // find cov adaptively
    /*
    BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                 GradientEstimator* gradient_estimator_in);
                                 */

    BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                 const std::vector<arma::mat> &covariances_in,
                                 GradientEstimator *gradient_estimator_in);

    BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                 const arma::mat &covariance_in,
                                 GradientEstimator *gradient_estimator_in);

    BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                 const arma::mat &covariance_in,
                                 double scale_in,
                                 GradientEstimator *gradient_estimator_in);

    BarkerDynamicsProposalKernel(const std::string &variable_name_in,
                                 const double &sd_in,
                                 GradientEstimator *gradient_estimator_in);

    /**
     * @brief Performs the barkerdynamicsproposalkernel operation.
     *
     * @param another The GradientEstimator instance to copy from.
     */
    BarkerDynamicsProposalKernel(const BarkerDynamicsProposalKernel &another);

    /**
     * @brief Assignment operator for GradientEstimator.
     *
     * @param another The GradientEstimator instance to copy from.
     */
    void operator=(const BarkerDynamicsProposalKernel &another);
    Kernel *duplicate() const;
    ProposalKernel *proposal_kernel_duplicate() const;

    /**
     * @brief Sets the proposal parameters.
     *
     * @param proposal_parameters_in The proposal parameters.
     */
    void set_proposal_parameters(Parameters *proposal_parameters_in);

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

    arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                                 const Particle &proposed_particle,
                                                 const Particle &old_particle);

    /*
     arma::mat specific_gradient_of_log(const std::string &variable,
     Particle &proposed_particle,
     Particle &old_particle,
     const Parameters &conditioned_on_parameters);
     */
    // virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
    //                                                      Particle &proposed_particle,
    //                                                      Particle &old_particle)=0;

    /*
     arma::mat specific_subsample_gradient_of_log(const std::string &variable,
     Particle &proposed_particle,
     Particle &old_particle,
     const Parameters &conditioned_on_parameters);
     */

    void make_copy(const BarkerDynamicsProposalKernel &another);

    boost::unordered_map<std::string, GaussianProposalInfo> proposal_info;

    boost::unordered_map<Particle *, ProposalStore> proposal_store;

    // must simulate each dimension of all variables specified in variable names independently from a symmetric proposal
    // to work, must work on transformed space, with transformation the same as this
    // stored here
    IndependentProposalKernel *proposal_simulate;

    // stored here
    GradientEstimator *gradient_estimator;

    // stored here
    Index *index;
  };
}

#endif
