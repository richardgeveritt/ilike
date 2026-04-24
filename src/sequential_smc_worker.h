#ifndef SEQUENTIALSMCWORKER_H
#define SEQUENTIALSMCWORKER_H

#include "smc_worker.h"

#include <vector>

#include "particles.h"

namespace ilike
{
  /**
   * @file sequential_smc_worker.h
   * @brief Defines the ParticleSimulator class.
   *
   * Provides particle simulator functionality.
   *
   * @namespace ilike
   * @class ParticleSimulator
   * @brief The particle simulator class.
   */


class ParticleSimulator;

class SequentialSMCWorker : public SMCWorker
{
public:
  
  /**
   * @brief Performs the sequentialsmcworker operation.
   */
  SequentialSMCWorker();
  /**
   * @brief Performs the ~sequentialsmcworker operation.
   */
  virtual ~SequentialSMCWorker();
  
  /**
   * @brief Performs the sequentialsmcworker operation.
   *
   * @param the_smc_in The the smc.
   */
  SequentialSMCWorker(SMC* the_smc_in);
  
  /**
   * @brief Performs the sequentialsmcworker operation.
   *
   * @param another The ParticleSimulator instance to copy from.
   */
  SequentialSMCWorker(const SequentialSMCWorker &another);
  /**
   * @brief Assignment operator for ParticleSimulator.
   *
   * @param another The ParticleSimulator instance to copy from.
   */
  void operator=(const SequentialSMCWorker &another);
  /**
   * @brief Creates a deep copy of this ParticleSimulator object.
   *
   * @return The result.
   */
  SMCWorker* duplicate() const;
  
  /**
   * @brief Performs the pf initial weight operation.
   *
   * @param current_particles The current particles.
   */
  void pf_initial_weight(Particles &current_particles);
  void weight(const Index* index,
              Particles &current_particles);
  void smcfixed_weight(const Index* index,
                       Particles &current_particles);
  void smcadaptive_given_smcfixed_weight(const Index* index,
                                         Particles &current_particles);
  void smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                  Particles &current_particles);
  void marginal_weight(const Index* index,
                       Particles &current_particles,
                       Particles &previous_particles,
                       ProposalKernel* proposal_kernel);
  void generic_weight(const Index* index,
                      Particles &current_particles,
                      Particles &previous_particles,
                      ProposalKernel* proposal_kernel,
                      ProposalKernel* L_kernel);
  void pf_weight(const Index* index,
                 Particles &current_particles,
                 Particles &previous_particles,
                 ProposalKernel* proposal_kernel);
  
  /*
   void weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void pf_initial_weight(Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void smcadaptive_given_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void smcadaptive_given_smcfixed_evaluate_target(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void marginal_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   void generic_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters);
   void pf_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_weight(const Index* index,
                        Particles &current_particles);
  /**
   * @brief Performs the subsample pf initial weight operation.
   *
   * @param current_particles The current particles.
   */
  void subsample_pf_initial_weight(Particles &current_particles);
  void subsample_smcfixed_weight(const Index* index,
                                 Particles &current_particles);
  void subsample_smcadaptive_given_smcfixed_weight(const Index* index,
                                                   Particles &current_particles);
  void subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
                                                            Particles &current_particles);
  void subsample_marginal_weight(const Index* index,
                                 Particles &current_particles,
                                 Particles &previous_particles,
                                 ProposalKernel* proposal_kernel);
  void subsample_generic_weight(const Index* index,
                                Particles &current_particles,
                                Particles &previous_particles,
                                ProposalKernel* proposal_kernel,
                                ProposalKernel* L_kernel);
  void subsample_pf_weight(const Index* index,
                           Particles &current_particles,
                           Particles &previous_particles,
                           ProposalKernel* proposal_kernel);
  
  /*
   void subsample_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void subsample_pf_initial_weight(Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void subsample_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void subsample_smcadaptive_given_smcfixed_weight(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void subsample_smcadaptive_given_smcfixed_evaluate_target(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters);
   void subsample_marginal_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   void subsample_generic_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   ProposalKernel* L_kernel,
   const Parameters &conditioned_on_parameters);
   void subsample_pf_weight(const Index* index,
   Particles &current_particles,
   Particles &previous_particles,
   ProposalKernel* proposal_kernel,
   const Parameters &conditioned_on_parameters);
   */
  
  //Particles simulated_particles() const;
  //Particles& simulated_particles();
  /**
   * @brief Returns the unnormalised log incremental weights.
   *
   * @return The result.
   */
  arma::colvec get_unnormalised_log_incremental_weights() const;
  
protected:
  
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param next_particles The next particles.
   */
  void specific_simulate(Particles* next_particles);
  //void specific_simulate_and_weight(const Parameters &conditioned_on_parameters);
  
  void specific_move(Particles* next_particles,
                     const Particles* current_particles);
  //void specific_simulate_and_weight(const Parameters &conditioned_on_parameters);
  
  
  void specific_simulate(Particles* next_particles,
                         const Parameters &conditioned_on_parameters);
  /*
   void specific_move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_specific_simulate(Particles* next_particles);
  
  void subsample_specific_move(Particles* next_particles,
                               const Particles* current_particles);
  
  void subsample_specific_simulate(Particles* next_particles,
                                   const Parameters &conditioned_on_parameters);
  /*
   void subsample_specific_move(Particles* next_particles,
   const Particles* current_particles,
   const Parameters &conditioned_on_parameters);
   */
  
  void make_copy(const SequentialSMCWorker &another);
  
  // Don't think we want to store this here, since found in particles (I think).
  //std::vector<Simulator*> simulate_priors;
  //std::vector<> simulate_for_likelihoods;
  
  // Not stored here.
  //Particles* particles;
  
  /** @brief The log unnormalised incremental weights. */
  arma::colvec log_unnormalised_incremental_weights;
  
};
}

#endif
