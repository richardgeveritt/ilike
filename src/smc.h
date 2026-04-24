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
  
  /** @brief The log likelihood file stream. */
  std::ofstream log_likelihood_file_stream;
  /** @brief The output lengths file stream. */
  std::ofstream output_lengths_file_stream;
  /** @brief The time file stream. */
  std::ofstream time_file_stream;
  /** @brief The vector variables file stream. */
  std::ofstream vector_variables_file_stream;
  /** @brief The vector variable sizes file stream. */
  std::ofstream vector_variable_sizes_file_stream;
  /** @brief The incremental log likelihood file stream. */
  std::ofstream incremental_log_likelihood_file_stream;
  /** @brief The resampled file stream. */
  std::ofstream resampled_file_stream;
  /** @brief The ess file stream. */
  std::ofstream ess_file_stream;
  /** @brief The schedule parameters file stream. */
  std::ofstream schedule_parameters_file_stream;
  /** @brief The vector points file stream. */
  std::ofstream vector_points_file_stream;
  /** @brief The ancestor index file stream. */
  std::ofstream ancestor_index_file_stream;
  /** @brief The any points file stream. */
  std::ofstream any_points_file_stream; // should be one for each member of Parameters
  /** @brief The normalised weights file stream. */
  std::ofstream normalised_weights_file_stream;
  /** @brief The unnormalised weights file stream. */
  std::ofstream unnormalised_weights_file_stream;
  
  //Parameters single_particle_is_step() const;
  
  //virtual void smc_step()=0;
  
  //virtual void weight_update()=0;
  
  //virtual double single_particle_weight_update() const=0;
  
  // Stored here.
  //SMCOutput* output;
  
};
}

#endif
