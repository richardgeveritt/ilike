#ifndef SMC_H
#define SMC_H

#include <functional>
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator.h"
#include "particles.h"
#include "distributions.h"
#include "parameters.h"
#include "sequencer.h"
#include "ilike_header.h"
#include "ilike_hdf5_utils.h"
#include <memory>

namespace ilike
{

class SMCOutput;
class SMCWorker;
class SequentialSMCWorker;
class SMCCriterion;
class MCMC;
class MoveOutput;
class ParticleSimulator;

class SMC : public LikelihoodEstimator
{
public:
  bool is_estimate() const { return true; } // (dev/plan.md 3.21 (c3))


  void set_write_factor_output(bool write_factor_output_in) { this->write_factor_output = write_factor_output_in; }
  // (dev/plan.md 3.21 (c5)) whether this SMC's output will be written; its
  // factors' estimators are told whether theirs will (written, and factor
  // outputs written)
  void set_outputs_written(bool written);
  void set_compact_stored_output(bool compact_stored_output_in) { this->compact_stored_output = compact_stored_output_in; }
  void set_chunk_size(size_t chunk_size_in) { this->chunk_size = chunk_size_in; }
  void set_force_rich_generations(bool force_rich_generations_in) { this->force_rich_generations = force_rich_generations_in; }
  void set_number_of_threads(size_t number_of_threads_in) { this->number_of_threads = number_of_threads_in; }
  void set_show_progress(bool show_progress_in, bool progress_counts_iterations_in)
  {
    this->show_progress = show_progress_in;
    this->progress_counts_iterations = progress_counts_iterations_in;
  }
  void set_waste_free_chains(size_t waste_free_chains_in) { this->waste_free_chains = waste_free_chains_in; }
  // (dev/plan.md 3.22 U4) the particles the next run starts from (as many as
  // the number of particles; used once). With from_first_target they are
  // draws from the first target, so the first weights are its ratio to the
  // next target (an annealed path from given states: exchange, SAV/MAV).
  void set_initial_particles(const std::vector<Parameters> &points, bool from_first_target);
  // (dev/plan.md 3.22 U4) an annealed path: parameters moved with the last
  // schedule variable, from start to end
  void set_path(const std::vector<std::string> &names, const Parameters &start, const Parameters &end)
  {
    this->sequencer.set_path(names, start, end);
  }
  // parallel = TRUE was asked for but some step ran sequentially (G4)
  bool parallel_was_refused() const { return this->parallel_refused; }
  // The index over which this SMC weights (for repointing a flat
  // generation's summaries to it, G5); NULL if none.
  virtual const Index* get_weighting_index() const { return NULL; }
  // True if, as an estimator inside a factor, this SMC reweights its single
  // generation at a target set by the caller each time (the importance
  // sampler of an ABC likelihood with a tolerance sequence, C62); then its
  // output's adaptive summary is that generation's summaries (G2/G5).
  virtual bool reweights_at_caller_target() const { return false; }
  // (dev/plan.md 3.21 (c2)) an inner SMC following the caller's target:
  // the random part of its step, at the caller's move stage
  virtual void advance_at_caller_target_smc(SMCOutput *current_state, const Parameters &conditioned_on_parameters) {}
  
  SMC();
  SMC(RandomNumberGenerator* rng_in,
      size_t* seed_in,
      Data* data_in,
      const Parameters &algorithm_parameters_in,
      size_t number_of_particles_in,
      size_t lag_in,
      size_t lag_proposed_in,
      const std::vector<const ProposalKernel*> proposals_in,
      //SummaryStatisticsPtr summary_statistics_in,
      double resampling_desired_ess_in,
      //EvaluateLogDistributionPtr evaluate_log_proposal_in,
      bool proposal_is_evaluated_in,
      bool smcfixed_flag_in,
      bool sequencer_limit_is_fixed_in,
      bool transform_proposed_particles,
      const std::string &results_name_in);
  
  SMC(RandomNumberGenerator* rng_in,
      size_t* seed_in,
      Data* data_in,
      const Parameters &algorithm_parameters_in,
      size_t number_of_particles_in,
      size_t lag_in,
      size_t lag_proposed_in,
      const std::vector<const ProposalKernel*> proposals_in,
      //SummaryStatisticsPtr summary_statistics_in,
      SMCCriterion* adaptive_resampling_method_in,
      //EvaluateLogDistributionPtr evaluate_log_proposal_in,
      bool proposal_is_evaluated_in,
      bool smcfixed_flag_in,
      bool sequencer_limit_is_fixed_in,
      bool transform_proposed_particles,
      const std::string &results_name_in);

  virtual ~SMC();
  
  SMC(const SMC &another);
  void operator=(const SMC &another);
  virtual SMC* smc_duplicate() const=0;
  
  SMCOutput* run();
  SMCOutput* run(const Parameters &conditioned_on_parameters);
  
  SMCOutput* run(const std::string &directory_name);
  SMCOutput* run(const std::string &directory_name,
                 const Parameters &conditioned_on_parameters);
  
  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  void setup();
  void setup(const Parameters &parameters);
  void setup_variables();
  void setup_variables(const Parameters &parameters);
  
  void resample(SMCOutput* current_state);
  
  virtual MoveOutput* move(RandomNumberGenerator &rng,
                           const Particle &particle) const=0;

  // Reuse of the objects of dropped generations (Stage 3 (a)): if true, the
  // output gives a dropped generation's objects to the next one, and the
  // worker calls move_reusing with one of them. move_reusing must return the
  // same as move(rng, particle), using spare (or deleting it).
  virtual bool reuses_particle_objects() const { return false; }
  virtual MoveOutput* move_reusing(RandomNumberGenerator &rng,
                                   const Particle &particle,
                                   MoveOutput* spare) const;

  virtual void weight_for_adapting_sequence(const Index* index,
                                            Particles &current_particles)=0;

  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     const Particle &particle) const=0;

  virtual void subsample_weight_for_adapting_sequence(const Index* index,
                                                      Particles &current_particles)=0;

protected:

  // Weighting at a trial target from adaptive summaries (dev/plan.md 2.6a,
  // G2): the incremental log weights of current_particles at the sequencer's
  // current schedule parameters, computed from the generation's summaries
  // when every factor of the index has one (otherwise the worker's rich
  // pass, as before). The weights are bit-for-bit those of the rich pass
  // (same expressions, same order of summation); only the particles' cached
  // factor values are not touched, so after the target has been chosen,
  // finish_weighting_from_summaries brings them to the state the rich passes
  // would have left them in (one rich pass at the last trial's values).
  // A generation simulated from the proposal in chunks and kept flat
  // (dev/plan.md 2.6a G1/G2, 2.6b G5; formerly ImportanceSampler's): each
  // point simulated, its smcfixed part over `index` evaluated, then either
  // its weight at the fixed target or its adaptive summary kept; its
  // variables (and factor records) kept if keep_points.
  void simulate_generation_in_chunks(SMCOutput *current_state,
                                     const Parameters* conditioned_on_parameters,
                                     bool summarise,
                                     const Index* index,
                                     bool keep_points);

  void weight_from_summaries_if_possible(const Index* index,
                                         Particles &current_particles);
  void finish_weighting_from_summaries(const Index* index,
                                       Particles &current_particles);
  
  void specific_change_data(Data* new_data);
  
  friend SequentialSMCWorker;
  
  bool proposal_is_evaluated;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<ParticleSimulator> particle_simulator;
  
  virtual SMCOutput* specific_run()=0;
  
  void simulate_proposal(SMCOutput* simulation);
  void finish_initial_particles(Particles &particles);
  
  virtual SMCOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  
  void simulate_proposal(SMCOutput* simulation,
                         const Parameters &conditioned_on_parameters);
  
  // Do what needs to be done for each particlar type of SMC (might iterate, might not).
  SMCOutput* initialise_smc();
  virtual SMCOutput* specific_initialise_smc()=0;
  virtual void simulate_smc(SMCOutput* simulation)=0;

  // (dev/plan.md 3.21) simulate as simulate_smc(simulation, parameters), but
  // stop as soon as the estimate is known to be below threshold (true: it
  // stopped, and the simulation is not usable); by default no early stop
  virtual bool simulate_smc_with_threshold(SMCOutput* simulation,
                                           const Parameters &conditioned_on_parameters,
                                           double threshold);
  // checked after each chunk of a sequential chunked generation: true stops it
  std::function<bool(Particles &generation, size_t points_done)> chunk_stop;
  // simulate and evaluate with the threshold (below: stopped early); by
  // default simulate_smc_with_threshold, then a full evaluation
  virtual double run_with_threshold(SMCOutput* simulation,
                                    const Parameters &conditioned_on_parameters,
                                    double threshold,
                                    bool &below);
  bool chunk_stopped = false;
  virtual void evaluate_smc(SMCOutput* simulation)=0;
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  SMCOutput* initialise_smc(const Parameters &parameters);
  virtual SMCOutput* specific_initialise_smc(const Parameters &parameters)=0;
  virtual void simulate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  virtual void evaluate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                          const Parameters &parameters)=0;
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                            const Parameters &parameters)=0;
  
  virtual void subsample_simulate_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  virtual void subsample_simulate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                                    const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                                      const Parameters &parameters)=0;
  
  void setup_default_ancestor_variables();

  void make_copy(const SMC &another);

  friend SMCWorker;
  friend SMCOutput;
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<SMCWorker> the_worker;
  // Stored here.
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<SMCCriterion> resampling_criterion;
  
  // Stored here.
  //Transform* transform;
  //bool store_raw;
  //bool store_transformed;
  
  Sequencer sequencer;
  // not stored here
  //Parameters* sequencer_parameters;
  bool sequencer_limit_is_fixed;
  size_t number_of_particles;
  size_t lag;
  size_t lag_proposed;

  bool proposed_particles_inputted;
  std::vector<Parameters> initial_particles; // not needed, unless initial values provided
  arma::colvec log_probabilities_of_initial_values;
  bool initial_particles_once = false;      // (U4: set_initial_particles, for the next run only)
  bool initial_from_first_target = false;
  
  std::vector<std::string> vector_variables;
  std::vector<std::string> any_variables;
  
  std::vector<size_t> vector_variable_sizes;
  
  std::vector<size_t> default_ancestor_variables;
  
  std::string results_name;
  
  bool initialised;
  
  // not stored here
  std::vector<const ProposalKernel*> proposals_to_transform_for;
  std::vector<const ProposalKernel*> proposals_to_find_gradient_for;
  
  std::shared_ptr<HighFive::File> h5_file;
  // per-step records of the algorithm (not of a generation), appended to the
  // output file when it is written (SMCOutput::write_to_file); none by default
public:
  virtual void write_step_records(HighFive::File &file) {}
protected:
  
  std::string h5_file_path;

  // Whether write_to_file writes each particle's factor outputs (the
  // iterationK/factorJ_* files). Set from the recipe (R argument
  // write_factor_output, default TRUE).
  bool write_factor_output = true;
  bool outputs_written = true;

  // Set for an SMC run used inside a factor (e.g. an algorithmic likelihood):
  // at the end of each run its stored generations are kept in compact form
  // and written concatenated over the outer points (dev/plan.md, nested
  // output design).
  bool compact_stored_output = false;

  // Chunked generations (dev/plan.md 2.6a): the number of points simulated
  // and weighted together with one set of reused particle objects and one
  // random-number stream (the R argument grain_size); 0 = the whole
  // generation. Results depend on it (not on whether a run is parallel).
  size_t chunk_size = 0;

  // Testing switch (recipe$force_rich_generations): never use the flat form.
  bool force_rich_generations = false;

  // Parallel runs (dev/plan.md 2.6a, G4): the R argument parallel, and the
  // number of threads to share the chunks between (the recipe's
  // number_of_threads, from RcppParallel's setting). Results are those of the
  // sequential run with the same seed and chunk size.
  bool parallel = false;
  size_t number_of_threads = 1;

  // Waste-free SMC (Dau & Chopin 2020; dev/plan.md 2.7): when > 0 (M), each
  // step resamples M ancestors and runs each as a chain of
  // number_of_particles/M states (the ancestor and that many minus one MCMC
  // moves), all of which form the next generation. 0: standard SMC.
  size_t waste_free_chains = 0;

  // Progress bar (dev/plan.md 3.11): only for the top-level algorithm (set
  // from the recipe by the R wrappers; never copied to inner estimators).
  // progress_counts_iterations: the bar counts MCMC iterations (MCMC()).
  // progress_counts_points: set while a chunked generation's points are what
  // the bar counts (IS).
  bool show_progress = false;
  bool progress_counts_iterations = false;
  bool progress_counts_points = false;
  bool parallel_refused = false;

  
};
}

#endif
