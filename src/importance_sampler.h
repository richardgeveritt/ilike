#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "smc.h"
#include "function_pointers.h"
//#include "model_and_algorithm.h"
#include "distributions.h"
#include "parameters.h"

class SMCOutput;
class MoveOutput;
class IndependentProposalKernel;

class ImportanceSampler : public SMC
{
public:

  ImportanceSampler();
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    SimulateIndependentProposalPtr simulate_prior_in,
                    bool smcfixed_flag_in,
                    bool parallel_in,
                    size_t grain_size_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                    SimulateIndependentProposalPtr simulate_proposal_in,
                    EvaluateLogDistributionPtr evaluate_log_proposal_in,
                    bool smcfixed_flag_in,
                    bool parallel_in,
                    size_t grain_size_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                    IndependentProposalKernel* proposal_in,
                    bool smcfixed_flag_in,
                    bool parallel_in,
                    size_t grain_size_in);

  ImportanceSampler(const ImportanceSampler &another);
  virtual ~ImportanceSampler(void);

  void operator=(const ImportanceSampler &another);
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
  
  SMCOutput* specific_run(const Parameters &conditioned_on_parameters);
  
  SMCOutput* initialise_smc(const Parameters &conditioned_on_parameters);

  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void subsample_simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                            const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                              const Parameters &conditioned_on_parameters);
  
  /*
  void weight(SMCOutput* simuation,
              const Parameters &conditioned_on_parameters);
  
  void smcfixed_weight(SMCOutput* simuation,
              const Parameters &conditioned_on_parameters);
  
  void smcadaptive_given_smcfixed_weight(SMCOutput* simuation,
                          const Parameters &conditioned_on_parameters);
  */
  //void smc_simulate(SMCOutput* current_state);
  //void smc_weight(SMCOutput* current_state);

  void make_copy(const ImportanceSampler &another);
  
  // stored here
  Index* index;

  // void smc_step();
  //
  // void weight_update();

  //double single_particle_weight_update() const;
};

#endif
