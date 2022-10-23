#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "smc.h"
#include "proposal_kernel.h"

class SMCOutput;
class MoveOutput;

class ParticleFilter : public SMC
{
public:
  
  // unsure

  ParticleFilter();

  ParticleFilter(RandomNumberGenerator* rng_in,
                 size_t* seed_in,
                 Data* data_in,
                 size_t number_of_particles_in,
                 size_t lag_in,
                 size_t lag_proposed_in,
                 double resampling_desired_ess_in,
                 ProposalKernel* proposal_kernel_in,
                 EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                 EvaluateLogDistributionPtr evaluate_log_prior_in,
                 SimulateDistributionPtr simulate_proposal_in,
                 EvaluateLogDistributionPtr evaluate_log_proposal_in,
                 bool parallel_in,
                 size_t grain_size_in);
  ParticleFilter(const ParticleFilter &another);
  virtual ~ParticleFilter(void);

  void operator=(const ParticleFilter &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
  void weight_for_adapting_sequence(Particles &current_particles);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle,
                   const Parameters &conditioned_on_parameters);
  
  void weight_for_adapting_sequence(Particles &current_particles,
                                    const Parameters &conditioned_on_parameters);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle,
                             const Parameters &conditioned_on_parameters);
  
  void subsample_weight_for_adapting_sequence(Particles &current_particles,
                                              const Parameters &conditioned_on_parameters);

protected:
  
  SMCOutput* specific_run();
  
  SMCOutput* initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  void evaluate_smc(SMCOutput* simulation);
  
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  //void mcmc_move(SMCOutput* current_state);
  
  
  SMCOutput* specific_run(const Parameters &parameters);
  
  SMCOutput* initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                    const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_simulate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                            const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                              const Parameters &conditioned_on_parameters);
  
  //void mcmc_move(SMCOutput* current_state,
  //               const Parameters &conditioned_on_parameters);

  //void smc_update(SMCOutput* current_state);

  void make_copy(const ParticleFilter &another);
  
  // stored here
  ProposalKernel* proposal_kernel;
  
  std::string index_name;
  //std::string time_name;
  std::vector<std::string> measurements_names;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;
  bool last_index_is_fixed;
  size_t lag;
  
  bool dynamic_proposal_is_evaluated;
  
  bool check_termination() const;

  //void smc_step();

  //void weight_update();
};

#endif
