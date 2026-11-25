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

  // (dev/plan.md 3.21) as a likelihood estimator in a Metropolis-Hastings
  // move: a bootstrap filter defers its run to evaluation, and stops when the
  // log-likelihood so far plus the measurement models' maxima at the
  // remaining steps is below the threshold
  bool defers_simulation() const;
  // (dev/plan.md 3.21 (c2'); Chopin, Jacob & Papaspiliopoulos's SMC^2) a
  // filter inside a factor that the caller's SMC takes in one observation at
  // a time (data tempering, data(fraction, …)): each output
  // keeps its particles between the caller's steps; an evaluation is a trial
  // (the current particles weighted by the next observation's density, no
  // random numbers); at the caller's move stage the output keeps that
  // weighting, resamples and propagates to the next time; a new output (a
  // new caller particle, or a PMMH proposal) is run through the observations
  // the caller has committed so far. With anneal (3.21 (c2''), F1) the
  // caller may take in a share r of the next observation: the trial is the
  // filter's estimate of the next observation's ratio to the power r (a
  // geometric bridge between the two extended targets), and the move stage
  // commits an observation only once it is whole
public:
  void set_follows_caller_target(const std::string &fraction_variable, bool anneal = false);
  bool advances_at_caller_target() const { return this->follows_caller_target; }
  size_t observations() const { return this->last_index - this->first_index + 1; }
  // (3.21 (c2''), F2) SMC^2's exchange step (Chopin, Jacob & Papaspiliopoulos
  // §3.6.1): the caller asks after its moves whether the acceptance rate calls
  // for more particles; if so it sets the new number and makes each of its
  // particles a new output (a filter of the new size through the same
  // observations), weighting it by the ratio of the two estimates
  void set_exchange(double below, double multiply, size_t maximum);
  bool has_exchange() const { return this->exchange_below >= 0.0; }
  size_t exchange_size(double acceptance_rate, const Parameters &caller_parameters) const; // (0: no exchange)
  void set_number_of_particles(size_t number_of_particles_in);
  size_t get_number_of_particles() const { return this->number_of_particles; }
  size_t get_number_of_inner_particles() const { return this->number_of_particles; }
  void set_number_of_inner_particles(size_t number_of_inner_particles_in) { this->set_number_of_particles(number_of_inner_particles_in); }
  void advance_at_caller_target_smc(SMCOutput *current_state, const Parameters &conditioned_on_parameters);
protected:
  bool follows_caller_target = false;
  std::string followed_variable;
  bool followed_anneal = false; // (F1: shares of an observation allowed)
  double exchange_below = -1.0; // (F2: no exchange when negative)
  double exchange_multiply = 2.0;
  size_t exchange_max = 0; // (no limit when 0)
  size_t followed_steps = 0; // (the most observations any output has committed)
  double initial_time = 0.0;
  void tempering_point(const Parameters &conditioned_on_parameters, size_t &k, double &r) const;
  double trial_value(SMCOutput *current_state, size_t k, double r);
  void set_time_step(size_t step);
  double weight_current(SMCOutput *current_state, double power = 1.0);
  void commit_observation(SMCOutput *current_state, const Parameters &conditioned_on_parameters);
  void follow_from_start(SMCOutput *current_state, const Parameters &conditioned_on_parameters);
  double run_with_threshold(SMCOutput* simulation,
                            const Parameters &conditioned_on_parameters,
                            double threshold,
                            bool &below);
  bool threshold_active = false;
  bool threshold_stopped = false;
  double stop_threshold = 0.0;
  std::vector<double> stop_remaining; // [k]: the sum of the bounds of the steps from first_index + k on

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
