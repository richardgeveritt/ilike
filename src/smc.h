#ifndef SMC_H
#define SMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include "model_and_algorithm.h"
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
  /**
   * @file smc.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


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

  void set_write_factor_output(bool write_factor_output_in) { this->write_factor_output = write_factor_output_in; }
  void set_compact_stored_output(bool compact_stored_output_in) { this->compact_stored_output = compact_stored_output_in; }
  void set_chunk_size(size_t chunk_size_in) { this->chunk_size = chunk_size_in; }
  void set_force_rich_generations(bool force_rich_generations_in) { this->force_rich_generations = force_rich_generations_in; }
  void set_number_of_threads(size_t number_of_threads_in) { this->number_of_threads = number_of_threads_in; }
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
  
  /**
   * @brief Performs the smc operation.
   */
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
  
  //SMC(RandomNumberGenerator* rng_in,
  //    size_t* seed_in,
  //    Data* data_in,
  //    size_t lag_in,
  //    size_t lag_proposed_in,
  //    double resampling_desired_ess_in,
  //    const std::vector<Parameters> &initial_values_in,
  //    const arma::colvec &log_probabilities_of_initial_values_in,
  //    bool smcfixed_flag_in,
  //    bool sequencer_limit_is_fixed_in);
  
  /**
   * @brief Performs the ~smc operation.
   */
  virtual ~SMC();
  
  /**
   * @brief Performs the smc operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  SMC(const SMC &another);
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const SMC &another);
  /**
   * @brief Creates a deep copy and returns it as a smc pointer.
   *
   * @return The result.
   */
  virtual SMC* smc_duplicate() const=0;
  
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @return The result.
   */
  SMCOutput* run();
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  SMCOutput* run(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @param directory_name The directory name.
   *
   * @return The result.
   */
  SMCOutput* run(const std::string &directory_name);
  SMCOutput* run(const std::string &directory_name,
                 const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  /**
   * @brief Performs the setup variables operation.
   */
  void setup_variables();
  /**
   * @brief Performs the setup variables operation.
   *
   * @param parameters The parameters.
   */
  void setup_variables(const Parameters &parameters);
  
  /**
   * @brief Resamples the particle population.
   *
   * @param current_state The current state.
   */
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
                                   MoveOutput* spare) const
  {
    delete spare;
    return this->move(rng, particle);
  }
  //virtual void specific_move(RandomNumberGenerator &rng,
  //                           Particle* new_particle,
  //                           Particle &particle)=0;
  
  //virtual void weight_for_adapting_sequence(const Index* index,
  //                                          Particles &current_particles)=0;
  
  /*
   virtual MoveOutput* move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters)=0;
   */
  //virtual void specific_move(RandomNumberGenerator &rng,
  //                           Particle* new_particle,
  //                           Particle &particle,
  //                           const Parameters &conditioned_on_parameters)=0;
  
  virtual void weight_for_adapting_sequence(const Index* index,
                                            Particles &current_particles)=0;
  
  /*
   virtual void weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
                                     const Particle &particle) const=0;
  
  /*
   virtual MoveOutput* subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  //virtual void specific_subsample_move(RandomNumberGenerator &rng,
  //                                     Particle* new_particle,
  //                                     Particle &particle,
  //                                     const Parameters &conditioned_on_parameters)=0;
  
  virtual void subsample_weight_for_adapting_sequence(const Index* index,
                                                      Particles &current_particles)=0;
  
  /*
   virtual void subsample_weight_for_adapting_sequence(const Index* index,
   Particles &current_particles,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  //EvaluateLogDistributionPtr evaluate_log_proposal;
  
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
  void weight_from_summaries_if_possible(const Index* index,
                                         Particles &current_particles);
  void finish_weighting_from_summaries(const Index* index,
                                       Particles &current_particles);
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  friend SequentialSMCWorker;
  friend class SimulateWorker;
  friend class ConditionalSimulateWorker;
  friend class SubsampleSimulateWorker;
  friend class SubsampleConditionalSimulateWorker;
  friend class MoveWorker;
  friend class SubsampleMoveWorker;
  
  /** @brief The proposal is evaluated. */
  bool proposal_is_evaluated;
  
  // stored here
  /** @brief The particle simulator. */
  ParticleSimulator* particle_simulator;
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  virtual SMCOutput* specific_run()=0;
  //virtual SMCOutput* specific_run(const std::string &directory_name)=0;
  
  /**
   * @brief Simulates proposal.
   *
   * @param simulation The simulation.
   */
  void simulate_proposal(SMCOutput* simulation);
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  virtual SMCOutput* specific_run(const Parameters &conditioned_on_parameters)=0;
  //virtual SMCOutput* specific_run(const std::string &directory_name,
  //                                const Parameters &conditioned_on_parameters)=0;
  
  void simulate_proposal(SMCOutput* simulation,
                         const Parameters &conditioned_on_parameters);
  
  // Do what needs to be done for each particlar type of SMC (might iterate, might not).
  /**
   * @brief Performs the initialise smc operation.
   *
   * @return The result.
   */
  SMCOutput* initialise_smc();
  /**
   * @brief Class-specific implementation for initialise smc.
   *
   * @return The result.
   */
  virtual SMCOutput* specific_initialise_smc()=0;
  /**
   * @brief Simulates smc.
   *
   * @param simulation The simulation.
   */
  virtual void simulate_smc(SMCOutput* simulation)=0;
  //virtual void evaluate_smc(SMCOutput* simulation,const Parameters &conditioned_on_parameters)=0;
  /**
   * @brief Evaluates the smc.
   *
   * @param simulation The simulation.
   */
  virtual void evaluate_smc(SMCOutput* simulation)=0;
  /**
   * @brief Evaluates the smcfixed part smc.
   *
   * @param simulation The simulation.
   */
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  /**
   * @brief Evaluates the smcadaptive part given smcfixed smc.
   *
   * @param simulation The simulation.
   */
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  /**
   * @brief Performs the initialise smc operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  SMCOutput* initialise_smc(const Parameters &parameters);
  /**
   * @brief Class-specific implementation for initialise smc.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  virtual SMCOutput* specific_initialise_smc(const Parameters &parameters)=0;
  virtual void simulate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  //virtual void evaluate_smc(SMCOutput* simulation,const Parameters &conditioned_on_parameters)=0;
  virtual void evaluate_smc(SMCOutput* simulation,
                            const Parameters &parameters)=0;
  virtual void evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                          const Parameters &parameters)=0;
  virtual void evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                            const Parameters &parameters)=0;
  
  /**
   * @brief Performs the subsample simulate smc operation.
   *
   * @param simulation The simulation.
   */
  virtual void subsample_simulate_smc(SMCOutput* simulation)=0;
  /**
   * @brief Performs the subsample evaluate smc operation.
   *
   * @param simulation The simulation.
   */
  virtual void subsample_evaluate_smc(SMCOutput* simulation)=0;
  /**
   * @brief Performs the subsample evaluate smcfixed part smc operation.
   *
   * @param simulation The simulation.
   */
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation)=0;
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed smc operation.
   *
   * @param simulation The simulation.
   */
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation)=0;
  
  virtual void subsample_simulate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smc(SMCOutput* simulation,
                                      const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcfixed_part_smc(SMCOutput* simulation,
                                                    const Parameters &parameters)=0;
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed_smc(SMCOutput* simulation,
                                                                      const Parameters &parameters)=0;
  
  /**
   * @brief Performs the setup default ancestor variables operation.
   */
  void setup_default_ancestor_variables();
  
  //virtual void smc_simulate(SMCOutput* current_state)=0;
  //virtual void smc_weight(SMCOutput* current_state)=0;
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const SMC &another);
  
  //Particles is_step() const;
  
  friend SMCWorker;
  friend SMCOutput;
  // Stored here.
  /** @brief The the worker. */
  SMCWorker* the_worker;
  // Stored here.
  /** @brief The resampling criterion. */
  SMCCriterion* resampling_criterion;
  
  // Stored here.
  //Transform* transform;
  //bool store_raw;
  //bool store_transformed;
  
  /** @brief The sequencer. */
  Sequencer sequencer;
  // not stored here
  //Parameters* sequencer_parameters;
  /** @brief The sequencer limit is fixed. */
  bool sequencer_limit_is_fixed;
  /** @brief The number of particles. */
  size_t number_of_particles;
  /** @brief The lag. */
  size_t lag;
  /** @brief The lag proposed. */
  size_t lag_proposed;
  
  //SummaryStatisticsPtr summary_statistics;
  
  /** @brief The proposed particles inputted. */
  bool proposed_particles_inputted;
  /** @brief The initial particles. */
  std::vector<Parameters> initial_particles; // not needed, unless initial values provided
  /** @brief The log probabilities of initial values. */
  arma::colvec log_probabilities_of_initial_values;
  
  /** @brief The vector variables. */
  std::vector<std::string> vector_variables;
  /** @brief The any variables. */
  std::vector<std::string> any_variables;
  
  /** @brief The vector variable sizes. */
  std::vector<size_t> vector_variable_sizes;
  
  /** @brief The default ancestor variables. */
  std::vector<size_t> default_ancestor_variables;
  
  /** @brief The results name. */
  std::string results_name;
  
  /** @brief The initialised. */
  bool initialised;
  
  // not stored here
  /** @brief The proposals to transform for. */
  std::vector<const ProposalKernel*> proposals_to_transform_for;
  /** @brief The proposals to find gradient for. */
  std::vector<const ProposalKernel*> proposals_to_find_gradient_for;
  
  /** @brief HDF5 output file (opened lazily in write_to_file). */
  std::shared_ptr<HighFive::File> h5_file;
  
  /** @brief Path to the HDF5 output file. */
  std::string h5_file_path;

  // Whether write_to_file writes each particle's factor outputs (the
  // iterationK/factorJ_* files). Set from the recipe (R argument
  // write_factor_output, default TRUE).
  bool write_factor_output = true;

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
  bool parallel_refused = false;

  //Parameters single_particle_is_step() const;
  
  //virtual void smc_step()=0;
  
  //virtual void weight_update()=0;
  
  //virtual double single_particle_weight_update() const=0;
  
  // Stored here.
  //SMCOutput* output;
  
};
}

#endif
