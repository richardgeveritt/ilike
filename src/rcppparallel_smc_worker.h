#ifndef RCPPPARALLELSMCWORKER_H
#define RCPPPARALLELSMCWORKER_H

#include <RcppParallel.h>

#include "smc_worker.h"

#include <vector>

#include "rcppparallel_workers.h"

namespace ilike
{
class SMC;

class RcppParallelSMCWorker : public SMCWorker
{
public:
  
  RcppParallelSMCWorker();
  RcppParallelSMCWorker(SMC* the_smc,
                        size_t grain_size_in);
  virtual ~RcppParallelSMCWorker();
  
  RcppParallelSMCWorker(const RcppParallelSMCWorker &another);
  void operator=(const RcppParallelSMCWorker &another);
  SMCWorker* duplicate() const;
  
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
  
  void subsample_weight(const Index* index,
                        Particles &current_particles);
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
  
  arma::colvec get_unnormalised_log_incremental_weights() const;
  
  //std::vector<Particle> get_particles() const;
  
  // Simulate from the proposal and weight.
  //void simulate_and_weight();
  
protected:
  
  // Just does the simulation from the proposal.
  void specific_simulate(Particles* next_particles);
  
  void specific_move(Particles* next_particles,
                     const Particles* current_particles);
  
  void specific_simulate(Particles* next_particles,
                         const Parameters &conditioned_on_parameters);
  
  void subsample_specific_simulate(Particles* next_particles);
  
  void subsample_specific_move(Particles* next_particles,
                               const Particles* current_particles);
  
  void subsample_specific_simulate(Particles* next_particles,
                                   const Parameters &conditioned_on_parameters);
  
  void make_copy(const RcppParallelSMCWorker &another);
  
  //uint64_t seed;
  //RandomNumberGenerator rng;
  
  size_t grain_size;
  
  friend SimulateWorker;
  friend ConditionalSimulateWorker;
  friend SubsampleSimulateWorker;
  friend SubsampleConditionalSimulateWorker;
  friend MoveWorker;
  friend SubsampleMoveWorker;
  
  SimulateWorker simulate_worker;
  ConditionalSimulateWorker conditional_simulate_worker;
  MoveWorker move_worker;
  WeightWorker weight_worker;
  PFInitialWeightWorker pf_initial_weight_worker;
  SMCFixedWeightWorker smcfixed_weight_worker;
  SMCAdaptiveGivenSMCFixedWeightWorker smcadaptive_given_smcfixed_weight_worker;
  SMCAdaptiveGivenSMCFixedEvaluateTargetWorker smcadaptive_given_smcfixed_evaluate_target_worker;
  MarginalWeightWorker marginal_weight_worker;
  GenericWeightWorker generic_weight_worker;
  PFWeightWorker pf_weight_worker;
  
  SubsampleSimulateWorker subsample_simulate_worker;
  SubsampleConditionalSimulateWorker subsample_conditional_simulate_worker;
  SubsampleMoveWorker subsample_move_worker;
  SubsampleWeightWorker subsample_weight_worker;
  SubsamplePFInitialWeightWorker subsample_pf_initial_weight_worker;
  SubsampleSMCFixedWeightWorker subsample_smcfixed_weight_worker;
  SubsampleSMCAdaptiveGivenSMCFixedWeightWorker subsample_smcadaptive_given_smcfixed_weight_worker;
  SubsampleSMCAdaptiveGivenSMCFixedEvaluateTargetWorker subsample_smcadaptive_given_smcfixed_evaluate_target_worker;
  SubsampleMarginalWeightWorker subsample_marginal_weight_worker;
  SubsampleGenericWeightWorker subsample_generic_weight_worker;
  SubsamplePFWeightWorker subsample_pf_weight_worker;
  
  std::vector<double>* log_unnormalised_incremental_weights;
  
};
}

#endif
