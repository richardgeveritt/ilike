#ifndef SMCMCMCMOVE_H
#define SMCMCMCMOVE_H

#include "particle_tuning.h"
#include "smc.h"
#include "ilike_header.h"
#include "independent_proposal_kernel.h"

namespace ilike
{

class ESJDAdaptor;

class SMCOutput;
class MoveOutput;
class Index;

class SMCMCMCMove : public SMC
{
public:
  
  SMCMCMCMove();
  
  // Multiple MCMC chains.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              const std::vector<Parameters> &initial_points_in,
              const arma::colvec &log_probabilities_of_initial_values_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // Multiple MCMC chains.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              IndependentProposalKernel* proposal_in,
              Index* without_cancelled_index,
              Index* full_index,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // Multiple MCMC chains.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              const std::vector<Parameters> &initial_points_in,
              const arma::colvec &log_probabilities_of_initial_values_in,
              Index* without_cancelled_index,
              Index* full_index,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // SMC with annealing.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  // SMC with annealing.
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              const std::vector<double> &temperatures_in,
              EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
              EvaluateLogDistributionPtr evaluate_log_prior_in,
              SimulateDistributionPtr simulate_proposal_in,
              EvaluateLogDistributionPtr evaluate_log_proposal_in,
              bool mcmc_at_last_step_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              double resampling_desired_ess_in,
              double annealing_desired_cess_in,
              size_t number_of_bisections_in,
              SMCTermination* termination_in,
              const std::string &sequence_variable_in,
              const std::vector<double> &schedule_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              IndependentProposalKernel* proposal_in,
              bool proposal_is_evaluated_in,
              bool smcfixed_flag_in,
              bool sequencer_limit_is_fixed_in,
              bool mcmc_at_last_step_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(RandomNumberGenerator* rng_in,
              size_t* seed_in,
              Data* data_in,
              const Parameters &algorithm_parameters,
              size_t number_of_particles_in,
              size_t lag_in,
              size_t lag_proposed_in,
              MCMC* mcmc_in,
              SMCCriterion* adaptive_resampling_in,
              SMCCriterion* adaptive_target_in,
              size_t number_of_bisections_in,
              SMCTermination* termination_in,
              const std::vector<std::string> &sequence_variable_in,
              const std::vector<std::vector<double>> &schedule_in,
              const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
              IndependentProposalKernel* proposal_in,
              Index* without_cancelled_index,
              Index* full_index,
              bool proposal_is_evaluated_in,
              bool smcfixed_flag_in,
              bool sequencer_limit_is_fixed_in,
              bool mcmc_at_last_step_in,
              bool transform_proposed_particles,
              bool parallel_in,
              size_t grain_size_in,
              const std::string &results_name_in);
  
  SMCMCMCMove(const SMCMCMCMove &another);
  
  virtual ~SMCMCMCMove();
  
  void operator=(const SMCMCMCMove &another);
  SMC* smc_duplicate() const;
  LikelihoodEstimator* duplicate() const;
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   const Particle &particle) const;

  void weight_for_adapting_sequence(const Index* index,
                                    Particles &current_particles);

  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             const Particle &particle) const;

  void subsample_weight_for_adapting_sequence(const Index* index,
                                              Particles &current_particles);

  // (dev/plan.md 3.21 (c2); Kerama, Thorne & Everitt, arXiv:2211.02172,
  // alg. 4) the inner SMC of integrate u with smc(…) inside a factor that the
  // caller's SMC sequence changes (rare-event ABC-SMC^2): its sequence
  // variable is the caller's (e.g. the ABC tolerance), and each output keeps
  // its points between the caller's steps. Evaluations are trials at the
  // caller's candidate targets (no random numbers); at the caller's move
  // stage it reweights to the chosen target, resamples and moves u
  // (advance_at_caller_target_smc); a new output (a new caller particle, or
  // a proposed one) is run through the caller's targets so far.
  void set_follows_caller_target(const std::string &variable);
  bool advances_at_caller_target() const { return this->follows_caller_target; }
  void advance_at_caller_target_smc(SMCOutput *current_state, const Parameters &conditioned_on_parameters);

protected:

  bool follows_caller_target = false;
  std::string followed_variable;
  std::vector<double> followed_targets; // the caller's targets so far
  double followed_target(const Parameters &conditioned_on_parameters) const;
  double step_to_target(SMCOutput *current_state, double target);
  void advance_at_current_target(SMCOutput *current_state, const Parameters &conditioned_on_parameters);
  void follow_from_start(SMCOutput *current_state, const Parameters &conditioned_on_parameters);
  void advance_inner_states(SMCOutput *current_state);

  SMCOutput* specific_run();

  SMCOutput* specific_initialise_smc();
  void simulate_smc(SMCOutput* simulation);
  void evaluate_smc(SMCOutput* simulation);
  void evaluate_smcfixed_part_smc(SMCOutput* simulation);
  void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation);
  
  void mcmc_move(SMCOutput* current_state);
  
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

  void make_copy(const SMCMCMCMove &another);

  // esjd tuning (dev/plan.md 3.19, C211): the kernels tuned (by instance
  // index) and their adaptors (inside this move's kernels), and each
  // particle's rows for the current move step (Particle::tuning points here;
  // N entries, kept in place, so the pointers stay valid)
  void prepare_tuning(SMCOutput* current_state);
public:
  // MCMC(…) runs its chains as one SMC move step: choose has no particles to
  // learn from there (C224)
  void set_chains_only(bool chains_only_in) { this->chains_only = chains_only_in; }
protected:
  bool chains_only = false;
  std::vector<std::pair<int, ESJDAdaptor*>> tuners;
  bool tuners_found = false;
  std::vector<ParticleTuning> particle_tunings;

  // MCMC moves on flat generations (dev/plan.md 2.6a, G3). When every factor
  // can be restored from its adaptive summary, the moved generation is kept
  // as numbers: every stored chain state (as the per-particle MCMC outputs
  // would have stored them), each chain's final state and that state's
  // adaptive summary. A move then loads its ancestor into one reused
  // particle object (from the numbers, or by copy from a rich generation),
  // runs the chain in place and appends the results. The random numbers and
  // every value are those of the per-particle path (one chunk); chunks of
  // chunk_size particles use their own streams, as in importance sampling.
  bool can_move_flat() const;
  bool can_simulate_flat() const;
  const Index* get_weighting_index() const { return this->index.get(); }

  // The per-particle (rich) move path in chunks (dev/plan.md 2.6b, G5): 0
  // not possible (the worker's path), 1 the factors use no random numbers
  // (shared by all chunks), 2 they include inner algorithms (one copy per
  // group of chunks, each particle's inner runs drawing from a seed counter
  // set from its chunk's stream).
  int rich_chunk_mode() const;
  void move_rich_in_chunks(Particles* next_particles,
                           Particles* current_particles,
                           bool copies);
  void move_flat(Particles* next_particles,
                 Particles* current_particles);
  // waste-free SMC (dev/plan.md 2.7): M ancestors, each run as a chain of
  // N/M states, all of them the next generation
  void move_waste_free(Particles* next_particles,
                       Particles* current_particles);
  void restore_particle(Particle &particle,
                        const Particles &generation,
                        size_t row,
                        const Parameters &restore_parameters) const;
  
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<MCMC> mcmc;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<Index> index;

  // Factors in the full target but not in index (they cancel in the SMC
  // weights, so the SMC never evaluates them). Evaluated before each MCMC run
  // so the moves' cached denominators are valid. Owned; may be NULL.
  std::unique_ptr<Index> cancelled_index;

  bool mcmc_at_last_step;

};
}

#endif
