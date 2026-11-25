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
  void set_unnormalised_log_incremental_weights(const arma::colvec &weights);
  
protected:
  
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
  
  void make_copy(const SequentialSMCWorker &another);
  
  // Don't think we want to store this here, since found in particles (I think).
  //std::vector<Simulator*> simulate_priors;
  //std::vector<> simulate_for_likelihoods;
  
  // Not stored here.
  //Particles* particles;
  
  arma::colvec log_unnormalised_incremental_weights;
public:
  // (dev/plan.md 3.21 (c2''), F2) after the SMC's number of particles changed
  void resize_for_particles() { this->log_unnormalised_incremental_weights.set_size(this->get_number_of_particles()); }
  
};
}

#endif
