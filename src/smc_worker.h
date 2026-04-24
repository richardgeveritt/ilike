#ifndef SMCWORKER_H
#define SMCWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "particles.h"
#include "distributions.h"
#include "ilike_header.h"
#include "proposal_kernel.h"
#include "sequencer.h"

namespace ilike
{
  /**
   * @file smc_worker.h
   * @brief Defines the SMC class.
   *
   * Provides smc functionality.
   *
   * @namespace ilike
   * @class SMC
   * @brief The smc class.
   */


class SMC;
class ParticleSimulator;
class LikelihoodEstimator;
class Index;

class SMCWorker
{
public:
  
  /**
   * @brief Performs the smcworker operation.
   */
  SMCWorker();
  
  /**
   * @brief Performs the smcworker operation.
   *
   * @param the_smc_in The the smc.
   */
  SMCWorker(SMC* the_smc_in);
  
  /**
   * @brief Performs the ~smcworker operation.
   *
   * @param void The void.
   */
  virtual ~SMCWorker(void);
  
  /**
   * @brief Performs the smcworker operation.
   *
   * @param another The SMC instance to copy from.
   */
  SMCWorker(const SMCWorker &another);
  /**
   * @brief Assignment operator for SMC.
   *
   * @param another The SMC instance to copy from.
   */
  void operator=(const SMCWorker &another);
  /**
   * @brief Creates a deep copy of this SMC object.
   *
   * @return The result.
   */
  virtual SMCWorker* duplicate() const=0;
  
  /**
   * @brief Returns the number of particles.
   *
   * @return The result.
   */
  size_t get_number_of_particles() const;
  /**
   * @brief Returns the rng.
   *
   * @return The result.
   */
  RandomNumberGenerator* get_rng();
  /**
   * @brief Returns the seed.
   *
   * @return The result.
   */
  size_t get_seed() const;
  /**
   * @brief Sets the seed.
   *
   * @param seed_in The seed.
   */
  void set_seed(size_t seed_in);
  
  /**
   * @brief Simulates the required variables.
   *
   * @param next_particles The next particles.
   */
  void simulate(Particles* next_particles);
  virtual void weight(const Index* index,
                      Particles &current_particles)=0;
  /**
   * @brief Performs the pf initial weight operation.
   *
   * @param current_particles The current particles.
   */
  virtual void pf_initial_weight(Particles &current_particles)=0;
  virtual void smcfixed_weight(const Index* index,
                               Particles &current_particles)=0;
  virtual void smcadaptive_given_smcfixed_weight(const Index* index,
                                                 Particles &current_particles)=0;
  virtual void smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                          Particles &current_particles)=0;
  virtual void marginal_weight(const Index* index,
                               Particles &current_particles,
                               Particles &previous_particles,
                               ProposalKernel* proposal_kernel)=0;
  virtual void generic_weight(const Index* index,
                              Particles &current_particles,
                              Particles &previous_particles,
                              ProposalKernel* proposal_kernel,
                              ProposalKernel* L_kernel)=0;
  virtual void pf_weight(const Index* index,
                         Particles &current_particles,
                         Particles &previous_particles,
                         ProposalKernel* proposal_kernel)=0;
  void move(Particles* next_particles,
            const Particles* current_particles);
  
  void simulate(Particles* next_particles,
                const Parameters &conditioned_on_parameters);
  
  /*
   virtual void weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void pf_initial_weight(Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void smcadaptive_given_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void smcadaptive_given_smcfixed_evaluate_target(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void marginal_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters)=0;
   virtual void generic_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters)=0;
   virtual void pf_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters)=0;
   void move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_simulate(Particles* next_particles);
  /**
   * @brief Performs the subsample pf initial weight operation.
   *
   * @param current_particles The current particles.
   */
  virtual void subsample_pf_initial_weight(Particles &current_particles)=0;
  virtual void subsample_weight(const Index* index,
                                Particles &current_particles)=0;
  virtual void subsample_smcfixed_weight(const Index* index,
                                         Particles &current_particles)=0;
  virtual void subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                           Particles &current_particles)=0;
  virtual void subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                                    Particles &current_particles)=0;
  virtual void subsample_marginal_weight(const Index* index,
                                         Particles &current_particles,
                                         Particles &previous_particles,
                                         ProposalKernel* proposal_kernel)=0;
  virtual void subsample_generic_weight(const Index* index,
                                        Particles &current_particles,
                                        Particles &previous_particles,
                                        ProposalKernel* proposal_kernel,
                                        ProposalKernel* L_kernel)=0;
  virtual void subsample_pf_weight(const Index* index,
                                   Particles &current_particles,
                                   Particles &previous_particles,
                                   ProposalKernel* proposal_kernel)=0;
  void subsample_move(Particles* next_particles,
                      const Particles* current_particles);
  
  void subsample_simulate(Particles* next_particles,
                          const Parameters &conditioned_on_parameters);
  /*
   virtual void subsample_pf_initial_weight(Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_smcadaptive_given_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_marginal_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_generic_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters)=0;
   virtual void subsample_pf_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters)=0;
   void subsample_move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  //void simulate_and_weight(const Parameters &conditioned_on_parameters);
  
  //virtual Particles& simulated_particles()=0;
  //virtual Particles simulated_particles() const=0;
  /**
   * @brief Returns the unnormalised log incremental weights.
   *
   * @return The result.
   */
  virtual arma::colvec get_unnormalised_log_incremental_weights() const=0;
  
protected:
  
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param next_particles The next particles.
   */
  virtual void specific_simulate(Particles* next_particles)=0;
  //virtual void specific_weight(const Parameters &conditioned_on_parameters)=0;
  //virtual void specific_simulate_and_weight(const Parameters &conditioned_on_parameters)=0;
  
  virtual void specific_move(Particles* next_particles,
                             const Particles* current_particles)=0;
  
  virtual void specific_simulate(Particles* next_particles,
                                 const Parameters &conditioned_on_parameters)=0;
  //virtual void specific_weight(const Parameters &conditioned_on_parameters)=0;
  //virtual void specific_simulate_and_weight(const Parameters &conditioned_on_parameters)=0;
  
  /*
   virtual void specific_move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual void subsample_specific_simulate(Particles* next_particles)=0;
  
  virtual void subsample_specific_move(Particles* next_particles,
                                       const Particles* current_particles)=0;
  
  virtual void subsample_specific_simulate(Particles* next_particles,
                                           const Parameters &conditioned_on_parameters)=0;
  
  /*
   virtual void subsample_specific_move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  void make_copy(const SMCWorker &another);
  
  friend Sequencer;
  /** @brief The the smc. */
  SMC* the_smc; // not stored here
  
  // stored in model_and_algorithm
  //Factors* factors;
  //ParticleSimulator* particle_simulator;
  //std::vector<LikelihoodEstimator*> likelihood_estimators;
  
  //Particles output;
  
  //EvaluateLogDistributionPtr evaluate_log_prior;
  
  //EvaluateLogDistributionPtr evaluate_log_proposal;
  
  // stored here
  // outer vector is which likelihood estimator
  // inner vector is which particle
  //std::vector< std::vector<LikelihoodEstimatorOutput*> > likelihood_estimator_outputs;
  
  // Add in function to do simulation of parts that cannot be done in parallel.
  
};
}

#endif
