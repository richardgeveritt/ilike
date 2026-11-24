#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "smc.h"
#include "ilike_header.h"
#include "distributions.h"
#include "parameters.h"

namespace ilike
{

  class SMCOutput;
  class MoveOutput;
  class IndependentProposalKernel;

  class ImportanceSampler : public SMC
  {
  public:
    ImportanceSampler();
    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                      SimulateDistributionPtr simulate_prior_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                      EvaluateLogDistributionPtr evaluate_log_prior_in,
                      SimulateDistributionPtr simulate_proposal_in,
                      EvaluateLogDistributionPtr evaluate_log_proposal_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      EvaluateLogLikelihoodPtr evaluate_log_likelihood_in,
                      EvaluateLogDistributionPtr evaluate_log_prior_in,
                      IndependentProposalKernel *proposal_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      const std::string &target_variable_in,
                      const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                      IndependentProposalKernel *proposal_in,
                      Index* without_cancelled_index,
                      bool proposal_is_evaluated_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool sequencer_limit_is_fixed_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      const std::string &target_variable_in,
                      double target_value_in,
                      const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                      IndependentProposalKernel *proposal_in,
                      Index* without_cancelled_index,
                      bool proposal_is_evaluated_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool sequencer_limit_is_fixed_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    // used in abc
    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      const std::string &target_variable_in,
                      const std::vector<std::string> &measurement_variables_in,
                      const std::string &scale_variable_in,
                      const arma::colvec &scale_in,
                      const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                      IndependentProposalKernel *proposal_in,
                      bool proposal_is_evaluated_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool sequencer_limit_is_fixed_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    // used in abc
    ImportanceSampler(RandomNumberGenerator *rng_in,
                      size_t *seed_in,
                      Data *data_in,
                      const Parameters &algorithm_parameters_in,
                      size_t number_of_particles_in,
                      const std::string &target_variable_in,
                      double target_value_in,
                      const std::vector<std::string> &measurement_variables_in,
                      const std::string &scale_variable_in,
                      const arma::colvec &scale_in,
                      const std::vector<LikelihoodEstimator *> &likelihood_estimators_in,
                      IndependentProposalKernel *proposal_in,
                      bool proposal_is_evaluated_in,
                      bool store_output,
                      bool smcfixed_flag_in,
                      bool sequencer_limit_is_fixed_in,
                      bool transform_proposed_particles,
                      bool parallel_in,
                      size_t grain_size_in,
                      const std::string &results_name_in);

    ImportanceSampler(const ImportanceSampler &another);
    virtual ~ImportanceSampler();

    void operator=(const ImportanceSampler &another);
    SMC *smc_duplicate() const;
    LikelihoodEstimator *duplicate() const;

    // G5 (dev/plan.md 2.6b): an inner importance sampler (e.g. an ABC
    // likelihood's) can be copied per thread when its generations are flat
    // (nothing of a run points to its rich particles) and not stored.
    bool supports_parallel_copies() const;
    const Index* get_weighting_index() const { return this->index.get(); }
    bool reweights_at_caller_target() const { return !this->sequencer_limit_is_fixed; }

    // (dev/plan.md 3.21, step (b)) inside a factor that the caller's SMC
    // sequence changes (integrate u with is(…) under an ABC tolerance
    // sequence): the points are kept, and each evaluation weights them at the
    // caller's current target (the factors read it through the shared
    // parameters), as in the ABC-SMC of arXiv:2211.02172 (alg. 3)
    void set_follows_caller_target(bool follows)
    {
      this->follows_caller_target = follows;
      this->smcfixed_flag = !follows;
    }

    MoveOutput *move(RandomNumberGenerator &rng,
                     const Particle &particle) const;

    void weight_for_adapting_sequence(const Index *index,
                                      Particles &current_particles);

    MoveOutput *subsample_move(RandomNumberGenerator &rng,
                               const Particle &particle) const;

    void subsample_weight_for_adapting_sequence(const Index *index,
                                                Particles &current_particles);

  protected:
    SMCOutput *specific_run();

    // Chunked importance sampling (dev/plan.md 2.6a, G1): when the target is
    // fixed, the points are simulated and weighted in chunks of chunk_size
    // with one set of reused particle objects, and the generation is kept in
    // flat form (numbers) instead of one rich particle per point.
    bool can_run_in_chunks() const;
    void simulate_and_weight_in_chunks(SMCOutput *current_state);

    // The same for a run conditioned on parameters whose target may change
    // after simulation (dev/plan.md 2.6a, G2: the importance sampler inside
    // an ABC likelihood, reweighted at every new tolerance): each point's
    // smcfixed part is evaluated as soon as it is simulated and its adaptive
    // summary kept, so that the generation is flat and every weighting is
    // done from the summaries.
    bool can_run_flat_conditioned() const;

    // (dev/plan.md 3.21) as a likelihood estimator with a single target and
    // a fixed scale (each point's weight final once it is simulated): its
    // run is deferred to evaluation, and may stop early
    bool defers_simulation() const;
    bool simulate_smc_with_threshold(SMCOutput* simulation,
                                     const Parameters &conditioned_on_parameters,
                                     double threshold);

    // The chunk loop. conditioned_on_parameters: NULL for an unconditioned
    // run. summarise: keep adaptive summaries (otherwise the weights at the
    // fixed target).
    void simulate_in_chunks(SMCOutput *current_state,
                            const Parameters* conditioned_on_parameters,
                            bool summarise);

    // Whether the chunks of a run may be processed in parallel (G4).
    bool can_run_in_parallel() const;

    SMCOutput *specific_initialise_smc();

    void simulate_smc(SMCOutput *simulation);
    void evaluate_smc(SMCOutput *simulation);
    void evaluate_smcfixed_part_smc(SMCOutput *simulation);
    void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *simulation);

    SMCOutput *specific_run(const Parameters &conditioned_on_parameters);

    SMCOutput *specific_initialise_smc(const Parameters &conditioned_on_parameters);

    void simulate_smc(SMCOutput *simulation,
                      const Parameters &conditioned_on_parameters);

    void evaluate_smc(SMCOutput *simulation,
                      const Parameters &conditioned_on_parameters);
    void evaluate_smcfixed_part_smc(SMCOutput *simulation,
                                    const Parameters &conditioned_on_parameters);
    void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *simulation,
                                                      const Parameters &conditioned_on_parameters);

    void subsample_simulate_smc(SMCOutput *simulation);

    void subsample_evaluate_smc(SMCOutput *simulation);
    void subsample_evaluate_smcfixed_part_smc(SMCOutput *simulation);
    void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *simulation);

    void subsample_simulate_smc(SMCOutput *simulation,
                                const Parameters &conditioned_on_parameters);
    void subsample_evaluate_smc(SMCOutput *simulation,
                                const Parameters &conditioned_on_parameters);
    void subsample_evaluate_smcfixed_part_smc(SMCOutput *simulation,
                                              const Parameters &conditioned_on_parameters);
    void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput *simulation,
                                                                const Parameters &conditioned_on_parameters);

    // This sets the scale of the abc kernel: function could be generalised to do other types of adaptation once the points are simulated
    void set_abc_scale(Particles* particles);
    void set_target_from_parameters(const Parameters &parameters);

    void make_copy(const ImportanceSampler &another);

    // the estimate at the caller's current target from the kept points: a
    // function of the points and the target only, so that a caller's
    // bisection over targets can try values in any order
    void weight_at_caller_target(SMCOutput *current_state, const Parameters &caller_parameters);
    bool follows_caller_target = false;

    // stored here
    // owned (dev/plan.md 3.5, CH3)
    std::unique_ptr<Index> index;

  };
}

#endif
