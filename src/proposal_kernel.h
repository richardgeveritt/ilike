#ifndef ProposalKernel_H
#define ProposalKernel_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "kernel.h"
#include "particle.h"
//#include "ensemble_member.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file proposal_kernel.h
   * @brief Defines the MCMCAdaptor class.
   *
   * An adaptor for mcmc tuning. Adjusts algorithm parameters based on accumulated statistics.
   *
   * @namespace ilike
   * @class MCMCAdaptor
   * @brief The mcmc adaptor class.
   */



class MCMCAdaptor;
class SMCAdaptor;
class Transform;

class SMCOutput;

class CompositeProposalKernel;
class CompositeIndependentProposalKernel;
class ReinforceGradientEstimator;
class DirectGradientEstimatorOutput;

class ProposalKernel : public Kernel
{
  
public:
  
  /**
   * @brief Performs the proposalkernel operation.
   */
  ProposalKernel();
  /**
   * @brief Performs the ~proposalkernel operation.
   */
  virtual ~ProposalKernel();
  
  /**
   * @brief Performs the proposalkernel operation.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  ProposalKernel(const Parameters &proposal_parameters_in);
  
  /**
   * @brief Performs the proposalkernel operation.
   *
   * @param another The MCMCAdaptor instance to copy from.
   */
  ProposalKernel(const ProposalKernel &another);
  
  /**
   * @brief Assignment operator for MCMCAdaptor.
   *
   * @param another The MCMCAdaptor instance to copy from.
   */
  void operator=(const ProposalKernel &another);
  /**
   * @brief Creates a deep copy and returns it as a proposal_kernel pointer.
   *
   * @return The result.
   */
  virtual ProposalKernel* proposal_kernel_duplicate() const=0;
  
  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  /*
   Particle move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  Particle subsample_move(RandomNumberGenerator &rng,
                          const Particle &particle) const;
  
  /*
   Particle subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  Particle subsample_move(RandomNumberGenerator &rng,
                          const std::string &variable,
                          const Particle &particle) const;
  
  /*
   Particle subsample_move(RandomNumberGenerator &rng,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  /*
   Particle move(RandomNumberGenerator &rng,
   const Index* index,
   Particle &particle) const;
   Particle move(RandomNumberGenerator &rng,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   //Particle subsample_move(RandomNumberGenerator &rng,
   //                              Particle &particle) const;
   Particle subsample_move(RandomNumberGenerator &rng,
   const Index* index,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   //Particle subsample_move(RandomNumberGenerator &rng,
   //                              const std::string &variable,
   //                              Particle &particle) const;
   Particle subsample_move(RandomNumberGenerator &rng,
   const Index* index,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
  /**
   * @brief Performs the smc adapt operation.
   *
   * @param current_state The current state.
   */
  void smc_adapt(SMCOutput* current_state);
  
  void mcmc_adapt(const Particle &current_particle,
                  size_t iteration_counter);
  
  //void use_transform(const Particle &particle);
  
  double evaluate_kernel(const Particle &proposed_particle,
                         const Particle &old_particle) const;
  double subsample_evaluate_kernel(const Particle &proposed_particle,
                                   const Particle &old_particle) const;
  
  //double evaluate_kernel_for_mh(const Particle &proposed_particle) const;
  //double subsample_evaluate_kernel_for_mh(const Particle &proposed_particle) const;
  
  /*
   double subsample_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const;
   */
  
  arma::mat gradient_of_log(const std::string &variable,
                            const Particle &proposed_particle,
                            const Particle &old_particle);
  
  /*
   arma::mat gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat subsample_gradient_of_log(const std::string &variable,
                                      const Particle &proposed_particle,
                                      const Particle &old_particle);
  
  /*
   arma::mat subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters);
   */
  
  virtual Parameters simulate(RandomNumberGenerator &rng,
                              const Particle &particle) const=0;
  
  /*
   virtual Parameters simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const Particle &particle) const=0;
  
  /*
   virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
                                        const std::string &variable,
                                        const Particle &particle) const=0;
  
  /*
   virtual Parameters subsample_simulate(RandomNumberGenerator &rng,
   const std::string &variable,
   Particle &particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  Transform* get_transform() const;
  
  /**
   * @brief Sets the proposal parameters.
   *
   * @param proposal_parameters_in The proposal parameters.
   */
  virtual void set_proposal_parameters(Parameters* proposal_parameters_in)=0;
  
  /**
   * @brief Simulates gradient estimator output.
   *
   * @return The result.
   */
  virtual GradientEstimatorOutput* simulate_gradient_estimator_output() const=0;
  
  /**
   * @brief Returns the proposals.
   *
   * @return The result.
   */
  virtual std::vector<const ProposalKernel*> get_proposals() const=0;
  
  /**
   * @brief Sets the index.
   *
   * @param index_in The index.
   */
  virtual void set_index(Index* index_in)=0;
  
  /**
   * @brief Sets the index if null.
   *
   * @param index_in The index.
   */
  virtual void set_index_if_null(Index* index_in)=0;
  
  /**
   * @brief Returns the instance index.
   *
   * @return The result.
   */
  int get_instance_index() const;
  
  /**
   * @brief Performs the can be evaluated operation.
   *
   * @return The result.
   */
  virtual bool can_be_evaluated() const=0;
  
  /**
   * @brief Sets the data.
   *
   * @param data_in The data.
   */
  virtual void set_data(Data* data_in)=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  friend CompositeProposalKernel;
  friend CompositeIndependentProposalKernel;
  friend ReinforceGradientEstimator;
  friend DirectGradientEstimatorOutput;
  
  virtual double specific_evaluate_kernel(const Particle &proposed_particle,
                                          const Particle &old_particle) const=0;
  
  /*
   virtual double specific_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                    const Particle &old_particle) const=0;
  
  /*
   virtual double specific_subsample_evaluate_kernel(Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters) const=0;
   */
  
  virtual arma::mat specific_gradient_of_log(const std::string &variable,
                                             const Particle &proposed_particle,
                                             const Particle &old_particle)=0;
  
  /*
   virtual arma::mat specific_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                                       const Particle &proposed_particle,
                                                       const Particle &old_particle)=0;
  
  /*
   virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
   Particle &proposed_particle,
   Particle &old_particle,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  void make_copy(const ProposalKernel &another);
  
  //EvaluateLogMCMCProposalPtr proposal_evaluate;
  
  //SimulateMCMCProposalPtr proposal_simulate;
  
  // stored here
  /** @brief The mcmc adaptor. */
  MCMCAdaptor* mcmc_adaptor;
  
  // stored here
  /** @brief The smc adaptor. */
  SMCAdaptor* smc_adaptor;
  
  // stored here
  /** @brief The transform. */
  Transform* transform;
  
  /** @brief The instance counter. */
  static int instance_counter; // Static member variable to track instances
  /** @brief The instance index. */
  int instance_index; // Instance-specific index
  
  //Parameters proposal_parameters;
  
};
}

#endif
