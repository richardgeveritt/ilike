#ifndef SEQUENTIALSMCWORKER_H
#define SEQUENTIALSMCWORKER_H

#include "smc_worker.h"

#include <vector>

#include "particles.h"

namespace ilike
{
class ParticleSimulator;

class SequentialSMCWorker : public SMCWorker
{
public:
  
  SequentialSMCWorker();
  virtual ~SequentialSMCWorker();
  
  SequentialSMCWorker(SMC* the_smc_in);
  
  SequentialSMCWorker(const SequentialSMCWorker &another);
  void operator=(const SequentialSMCWorker &another);
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
  arma::colvec get_unnormalised_log_incremental_weights() const;
  
protected:
  
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
  
  arma::colvec log_unnormalised_incremental_weights;
  
};
}

#endif
