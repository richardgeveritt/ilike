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

class SMC;
class ParticleSimulator;
class LikelihoodEstimator;
class Index;

class SMCWorker
{
public:
  
  SMCWorker();
  
  SMCWorker(SMC* the_smc_in);

  // The SMC this worker works for (after the SMC has been copied, B29).
  void set_smc(SMC* the_smc_in) { this->the_smc = the_smc_in; }
  
  virtual ~SMCWorker(void);
  
  SMCWorker(const SMCWorker &another);
  void operator=(const SMCWorker &another);
  virtual SMCWorker* duplicate() const=0;
  
  size_t get_number_of_particles() const;
  RandomNumberGenerator* get_rng();
  size_t get_seed() const;
  void set_seed(size_t seed_in);

  // Use these as the incremental log weights of the current weighting (for a
  // generation whose weights were computed while it was simulated in chunks,
  // dev/plan.md 2.6a). Default: not supported.
  virtual void set_unnormalised_log_incremental_weights(const arma::colvec &weights)
  {
    Rcpp::stop("SMCWorker::set_unnormalised_log_incremental_weights: not available for this worker.");
  }
  
  void simulate(Particles* next_particles);
  virtual void weight(const Index* index,
                      Particles &current_particles)=0;
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

  void subsample_simulate(Particles* next_particles);
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
  
  virtual arma::colvec get_unnormalised_log_incremental_weights() const=0;
  
protected:
  
  virtual void specific_simulate(Particles* next_particles)=0;
  
  virtual void specific_move(Particles* next_particles,
                             const Particles* current_particles)=0;
  
  virtual void specific_simulate(Particles* next_particles,
                                 const Parameters &conditioned_on_parameters)=0;

  virtual void subsample_specific_simulate(Particles* next_particles)=0;
  
  virtual void subsample_specific_move(Particles* next_particles,
                                       const Particles* current_particles)=0;
  
  virtual void subsample_specific_simulate(Particles* next_particles,
                                           const Parameters &conditioned_on_parameters)=0;

  void make_copy(const SMCWorker &another);
  
  friend Sequencer;
  SMC* the_smc; // not stored here
  
  // stored in model_and_algorithm
  //Factors* factors;
  //ParticleSimulator* particle_simulator;
  //std::vector<LikelihoodEstimator*> likelihood_estimators;

  // stored here
  // outer vector is which likelihood estimator
  // inner vector is which particle
  //std::vector< std::vector<LikelihoodEstimatorOutput*> > likelihood_estimator_outputs;
  
  // Add in function to do simulation of parts that cannot be done in parallel.
  
};
}

#endif
