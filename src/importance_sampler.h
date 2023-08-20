#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "smc.h"
#include "ilike_header.h"
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
                    const Parameters &algorithm_parameters_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    SimulateDistributionPtr simulate_prior_in,
                    bool smcfixed_flag_in,
                    bool transform_proposed_particles,
                    bool parallel_in,
                    size_t grain_size_in,
                    const std::string &results_name_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    const Parameters &algorithm_parameters_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                    SimulateDistributionPtr simulate_proposal_in,
                    EvaluateLogDistributionPtr evaluate_log_proposal_in,
                    bool smcfixed_flag_in,
                    bool transform_proposed_particles,
                    bool parallel_in,
                    size_t grain_size_in,
                    const std::string &results_name_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    const Parameters &algorithm_parameters_in,
                    size_t number_of_particles_in,
                    EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                    EvaluateLogDistributionPtr evaluate_log_prior_in,
                    IndependentProposalKernel* proposal_in,
                    bool smcfixed_flag_in,
                    bool transform_proposed_particles,
                    bool parallel_in,
                    size_t grain_size_in,
                    const std::string &results_name_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    const Parameters &algorithm_parameters_in,
                    size_t number_of_particles_in,
                    const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                    IndependentProposalKernel* proposal_in,
                    bool proposal_is_evaluated_in,
                    bool smcfixed_flag_in,
                    bool sequencer_limit_is_fixed_in,
                    bool transform_proposed_particles,
                    bool parallel_in,
                    size_t grain_size_in,
                    const std::string &results_name_in);
  
  ImportanceSampler(RandomNumberGenerator* rng_in,
                    size_t* seed_in,
                    Data* data_in,
                    const Parameters &algorithm_parameters_in,
                    size_t number_of_particles_in,
                    const std::string &target_variable_in,
                    double target_value_in,
                    const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                    IndependentProposalKernel* proposal_in,
                    bool proposal_is_evaluated_in,
                    bool smcfixed_flag_in,
                    bool sequencer_limit_is_fixed_in,
                    bool transform_proposed_particles,
                    bool parallel_in,
                    size_t grain_size_in,
                    const std::string &results_name_in);

  ImportanceSampler(const ImportanceSampler &another);
  virtual ~ImportanceSampler();

  void operator=(const ImportanceSampler &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle);
  
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
  
  SMCOutput* specific_run(const Parameters &conditioned_on_parameters);
  //SMCOutput* specific_run(const std::string &directory_name,
  //                        const Parameters &conditioned_on_parameters);
  
  SMCOutput* specific_initialise_smc(const Parameters &conditioned_on_parameters);

  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void subsample_simulate_smc(SMCOutput* simulation);
  
  void subsample_evaluate_smc(SMCOutput* simulation);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
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
