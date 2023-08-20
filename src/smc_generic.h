#ifndef SMCGENERIC_H
#define SMCGENERIC_H

#include "smc.h"
#include "proposal_kernel.h"

class SMCOutput;
class MoveOutput;
class Index;

class SMCGeneric : public SMC
{
public:

  SMCGeneric();

  SMCGeneric(RandomNumberGenerator* rng_in,
             size_t* seed_in,
             Data* data_in,
             size_t number_of_particles_in,
             size_t lag_in,
             size_t lag_proposed_in,
             double resampling_desired_ess_in,
             const std::vector<double> &temperatures_in,
             ProposalKernel* proposal_kernel_in,
             ProposalKernel* L_kernel_in,
             EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
             EvaluateLogDistributionPtr evaluate_log_prior_in,
             SimulateDistributionPtr simulate_proposal_in,
             EvaluateLogDistributionPtr evaluate_log_proposal_in,
             bool transform_proposed_particles,
             bool parallel_in,
             size_t grain_size_in,
             const std::string &results_name_in);
  SMCGeneric(const SMCGeneric &another);
  virtual ~SMCGeneric();

  void operator=(const SMCGeneric &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle);
  
  //void weight_for_adapting_sequence(Particles &current_particles);
  
  /*
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle,
                   const Parameters &conditioned_on_parameters);
  */
  
  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles);
  
  /*
  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles,
                                    const Parameters &conditioned_on_parameters);
  */
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             const Particle &particle);
  
  /*
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle,
                             const Parameters &conditioned_on_parameters);
  */
  
  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles);
  
  /*
  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles,
                                              const Parameters &conditioned_on_parameters);
  */

protected:
  
  SMCOutput* specific_run();
  //SMCOutput* specific_run(const std::string &directory_name);
  
  SMCOutput* specific_initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  
  void evaluate_smc(SMCOutput* simulation);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void mcmc_move(SMCOutput* current_state);
  
  SMCOutput* specific_run(const Parameters &parameters);
  /*
  SMCOutput* specific_run(const std::string &directory_name,
                          const Parameters &parameters);
  */
  SMCOutput* specific_initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  //void mcmc_move(SMCOutput* current_state,
  //               const Parameters &conditioned_on_parameters);
  
  void subsample_simulate_smc(SMCOutput* simulation);
  
  void subsample_simulate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate_smc(SMCOutput* simulation);
  void subsample_evaluate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                            const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                              const Parameters &conditioned_on_parameters);

  //void smc_update(SMCOutput* current_state);

  void make_copy(const SMCGeneric &another);
  
  // stored here
  ProposalKernel* proposal_kernel;
  ProposalKernel* L_kernel;
  
  // Stored here.
  //SMCCriterion* smc_criterion;
  
  // stored here
  Index* index;

  //void smc_step();

  //void weight_update();
};

#endif
