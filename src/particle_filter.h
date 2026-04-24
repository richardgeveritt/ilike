#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "smc.h"

namespace ilike
{
  /**
   * @file particle_filter.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class MoveOutput;
class IndependentProposalKernel;
class HMMIndex;
class VectorIndex;

class ParticleFilter : public SMC
{
public:
  
  // unsure
  
  /**
   * @brief Performs the particlefilter operation.
   */
  ParticleFilter();
  
  ParticleFilter(RandomNumberGenerator* rng_in,
                 size_t* seed_in,
                 Data* data_in,
                 const Parameters &algorithm_parameters,
                 size_t number_of_particles_in,
                 size_t lag_in,
                 size_t lag_proposed_in,
                 const std::string &index_name_in,
                 const std::string &time_name_in,
                 const std::string &time_diff_name_in,
                 size_t first_index_in,
                 size_t last_index_in,
                 size_t predictions_per_update_in,
                 double update_time_step_in,
                 double initial_time_in,
                 SMCCriterion* adaptive_resampling_in,
                 const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                 IndependentProposalKernel* proposal_in,
                 ProposalKernel* transition_model_in,
                 ProposalKernel* transition_proposal_in,
                 VectorIndex* evaluated_in_initial_weight_update,
                 VectorIndex* evaluated_in_pf_weight_update,
                 bool proposal_is_evaluated_in,
                 bool transition_proposal_is_evaluated_in,
                 bool smcfixed_flag_in,
                 bool sequencer_limit_is_fixed_in,
                 bool transform_proposed_particles,
                 bool parallel_in,
                 size_t grain_size_in,
                 const std::string &results_name_in);
  /**
   * @brief Performs the particlefilter operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  ParticleFilter(const ParticleFilter &another);
  /**
   * @brief Performs the ~particlefilter operation.
   *
   * @param void The void.
   */
  virtual ~ParticleFilter(void);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const ParticleFilter &another);
  /**
   * @brief Creates a deep copy and returns it as a smc pointer.
   *
   * @return The result.
   */
  SMC* smc_duplicate() const;
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle) const;
  
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
                             const Particle &particle) const;
  
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
  
  /**
   * @brief Performs the increment time index operation.
   */
  void increment_time_index();
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  SMCOutput* specific_run();
  //SMCOutput* specific_run(const std::string &directory_name);
  
  /**
   * @brief Class-specific implementation for initialise smc.
   *
   * @return The result.
   */
  SMCOutput* specific_initialise_smc();
  /**
   * @brief Simulates smc.
   *
   * @param simulation The simulation.
   */
  void simulate_smc(SMCOutput* simulation);
  /**
   * @brief Evaluates the smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smc(SMCOutput* simulation);
  
  /**
   * @brief Evaluates the smcfixed part smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  /**
   * @brief Evaluates the smcadaptive part given smcfixed smc.
   *
   * @param simulation The simulation.
   */
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  //void mcmc_move(SMCOutput* current_state);
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  SMCOutput* specific_run(const Parameters &parameters);
  //SMCOutput* specific_run(const std::string &directory_name,
  //                        const Parameters &parameters);
  
  /**
   * @brief Class-specific implementation for initialise smc.
   *
   * @param parameters The parameters.
   *
   * @return The result.
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
  
  /**
   * @brief Performs the subsample evaluate smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smc(SMCOutput* simulation);
  /**
   * @brief Performs the subsample simulate smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_simulate_smc(SMCOutput* simulation);
  /**
   * @brief Performs the subsample evaluate smcfixed part smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed smc operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
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
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const ParticleFilter &another);
  
  // stored here
  /** @brief The transition proposal. */
  ProposalKernel* transition_proposal;
  
  
  /** @brief The index name. */
  std::string index_name;
  /** @brief The time name. */
  std::string time_name;
  /** @brief The time diff name. */
  std::string time_diff_name;
  //std::vector<std::string> measurements_names;
  /** @brief The first index. */
  size_t first_index;
  /** @brief The last index. */
  size_t last_index;
  /** @brief The predictions per update. */
  size_t predictions_per_update;
  /** @brief The update time step. */
  double update_time_step;
  /** @brief The current time. */
  double current_time;
  /** @brief The current index. */
  size_t current_index;
  
  //std::vector<std::string> measurements_names;
  //size_t lag;
  
  // stored here
  /** @brief The index. */
  HMMIndex* index;
  
  /** @brief The transition proposal is evaluated. */
  bool transition_proposal_is_evaluated;
  
  /**
   * @brief Performs the check termination operation.
   *
   * @return The result.
   */
  bool check_termination() const;
  
  //void smc_step();
  
  //void weight_update();
};
}

#endif
