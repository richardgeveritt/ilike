#ifndef PARTICLEFILTER_H
#define PARTICLEFILTER_H

#include "smc.h"

namespace ilike
{

class SMCOutput;
class MoveOutput;
class IndependentProposalKernel;
class HMMIndex;
class VectorIndex;

class ParticleFilter : public SMC
{
public:
  
  // unsure
  
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
  ParticleFilter(const ParticleFilter &another);
  virtual ~ParticleFilter(void);
  
  void operator=(const ParticleFilter &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;

  // G5 (dev/plan.md 2.6b): an inner particle filter can be copied per
  // thread when what it keeps of a run outlives no rich particle (its stored
  // generations are compacted, or none is stored) and its data model can
  // be copied too.
  bool supports_parallel_copies() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle) const;

  bool reuses_particle_objects() const;
  
  MoveOutput* move_reusing(RandomNumberGenerator &rng,
  
                           const Particle &particle,
  
                           MoveOutput* spare) const;

  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles);

  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             const Particle &particle) const;

  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles);

protected:
  
  void increment_time_index();

  // Chunked moves and weights (dev/plan.md 2.6b, G5): when the data model
  // uses no random numbers and keeps no run state, the particles are moved
  // grain_size at a time, each chunk with its own random-number stream, and
  // weighted by chunks; in parallel when parallel = TRUE. The transition
  // kernels are only read (their data is set once per step, before the
  // chunks).
  bool can_move_in_chunks() const;
  void move_in_chunks(Particles* next_particles,
                      Particles* current_particles);
  void weight_in_chunks(SMCOutput* current_state);
  MoveOutput* move_with_data_set(RandomNumberGenerator &rng,
                                 const Particle &particle) const;
  MoveOutput* move_reusing_with_data_set(RandomNumberGenerator &rng,
                                         const Particle &particle,
                                         MoveOutput* spare) const;
  
  SMCOutput* specific_run();
  
  SMCOutput* specific_initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  void evaluate_smc(SMCOutput* simulation);
  
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  SMCOutput* specific_run(const Parameters &parameters);
  
  SMCOutput* specific_initialise_smc(const Parameters &parameters);
  void simulate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smc(SMCOutput* simulation,
                    const Parameters &conditioned_on_parameters);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                  const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                    const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate_smc(SMCOutput* simulation);
  void subsample_simulate_smc(SMCOutput* simulation);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void subsample_evaluate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_simulate_smc(SMCOutput* simulation,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                            const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                              const Parameters &conditioned_on_parameters);

  void make_copy(const ParticleFilter &another);
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<ProposalKernel> transition_proposal;

  std::string index_name;
  std::string time_name;
  std::string time_diff_name;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;

  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<HMMIndex> index;
  
  bool transition_proposal_is_evaluated;
  
  bool check_termination() const;

};
}

#endif
