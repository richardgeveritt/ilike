#ifndef ITERATIVEENSEMBLEKALMANINVERSION_H
#define ITERATIVEENSEMBLEKALMANINVERSION_H

#include "ensemble_kalman.h"
#include "ensemble_sequencer.h"

namespace ilike
{
  /**
   * @file ensemble_kalman_inversion.h
   * @brief Defines the EnsembleKalmanOutput class.
   *
   * Stores and manages the output produced by EnsembleKalman. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class EnsembleKalmanOutput
   * @brief The ensemble kalman output class.
   */


class EnsembleKalmanOutput;
class MCMC;
class MoveOutput;
class GaussianNoiseProposalKernel;

class EnsembleKalmanInversion : public EnsembleKalman
{
public:
  
  /**
   * @brief Performs the ensemblekalmaninversion operation.
   */
  EnsembleKalmanInversion();
  
  // single move from start to end
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          IndependentProposalKernel* prior_in,
                          SimulateModelPtr simulate_model_in,
                          std::shared_ptr<Transform> summary_statistics_in,
                          std::shared_ptr<Transform> transform_in,
                          const std::vector<std::string> &measurement_variables_in,
                          double significance_level_in,
                          size_t estimator_type_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  // adaptive annealing
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          double annealing_desired_cess_in,
                          size_t number_of_bisections_in,
                          const std::string &sequence_variable_in,
                          const std::vector<double> &schedule_in,
                          IndependentProposalKernel* prior_in,
                          SimulateModelPtr simulator_model_in,
                          std::shared_ptr<Transform> summary_statistics_in,
                          std::shared_ptr<Transform> transform_in,
                          const std::vector<std::string> &measurement_variables_in,
                          double significance_level_in,
                          size_t estimator_type_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  // adaptive annealing
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          double annealing_desired_cess_in,
                          size_t number_of_bisections_in,
                          const std::string &sequence_variable_in,
                          const std::vector<double> &schedule_in,
                          IndependentProposalKernel* prior_in,
                          std::shared_ptr<Transform> measurement_transform_function_in,
                          const std::vector<std::string> &measurement_variables,
                          const std::vector<GetMatrixPtr> &measurement_noise_functions,
                          std::shared_ptr<Transform> summary_statistics_in,
                          std::shared_ptr<Transform> transform_in,
                          double significance_level_in,
                          size_t estimator_type_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  // abc
  /*
   EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
   size_t* seed_in,
   Data* data_in,
   size_t number_of_ensemble_members_in,
   size_t lag_in,
   EnsembleShifter* shifter_in,
   double annealing_desired_cess_in,
   size_t number_of_bisections_in,
   const std::string &sequence_variable_in,
   const std::vector<double> &schedule_in,
   IndependentProposalKernel* prior_in,
   const std::vector<std::string> &measurement_variables_in,
   double min_epsilon_in,
   const std::string &scale_variable_in,
   std::shared_ptr<Transform> summary_statistics_in,
   std::shared_ptr<Transform> transform_in,
   double significance_level_in,
   bool parallel_in,
   size_t grain_size_in,
   const std::string &results_name_in);
   */
  
  // abc
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          size_t number_of_targets,
                          //size_t number_of_bisections_in,
                          const std::string &sequence_variable_in,
                          IndependentProposalKernel* prior_in,
                          const std::vector<std::string> &measurement_variables_in,
                          double min_epsilon_in,
                          const std::string &scale_variable_in,
                          const arma::colvec &scale_in,
                          std::shared_ptr<Transform> summary_statistics_in,
                          std::shared_ptr<Transform> transform_in,
                          double significance_level_in,
                          size_t estimator_type_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          SMCCriterion* adaptive_target_criterion,
                          size_t number_of_bisections_in,
                          SMCTermination* termination_in,
                          const std::string &sequence_variable_in,
                          const std::vector<double> &schedule_in,
                          IndependentProposalKernel* prior_in,
                          const std::vector<LikelihoodEstimator*> &likelihood_estimators_in,
                          const std::vector<MeasurementCovarianceEstimator*> &estimators_in,
                          std::shared_ptr<Transform> transform_in,
                          double significance_level_in,
                          size_t estimator_type_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  /**
   * @brief Performs the ensemblekalmaninversion operation.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  EnsembleKalmanInversion(const EnsembleKalmanInversion &another);
  
  /**
   * @brief Performs the ~ensemblekalmaninversion operation.
   */
  virtual ~EnsembleKalmanInversion();
  
  /**
   * @brief Assignment operator for EnsembleKalmanOutput.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void operator=(const EnsembleKalmanInversion &another);
  /**
   * @brief Creates a deep copy and returns it as a ensemble_kalman pointer.
   *
   * @return The result.
   */
  EnsembleKalman* ensemble_kalman_duplicate() const;
  /**
   * @brief Creates a deep copy of this EnsembleKalmanOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  /**
   * @brief Sets the abc schedule.
   *
   * @param current_state The current state.
   */
  void set_abc_schedule(EnsembleKalmanOutput* current_state);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle);
  
  /*
   MoveOutput* move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   
   MoveOutput* subsample_move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  //void weight_for_adapting_sequence(Particles &current_particles);
  
  //void weight_for_adapting_sequence(Particles &current_particles,
  //                                  const Parameters &conditioned_on_parameters);
  
  //void subsample_weight_for_adapting_sequence(Particles &current_particles,
  //                                            const Parameters &conditioned_on_parameters);
  
  //void weight_for_adapting_sequence(Ensemble &current_particles,
  //                                  double incremental_temperature);
  
protected:
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_run();
  
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise();
  /**
   * @brief Performs the ensemble kalman simulate operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation);
  /**
   * @brief Performs the ensemble kalman evaluate operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation);
  /**
   * @brief Performs the ensemble kalman evaluate smcfixed part operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation);
  /**
   * @brief Performs the ensemble kalman evaluate smcadaptive part given smcfixed operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation);
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_run(const Parameters &parameters);
  
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters);
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the ensemble kalman subsample simulate operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                          const Parameters &conditioned_on_parameters);
  
  //void smc_update(EnsembleKalmanOutput* current_state);
  
  //void setup_variables();
  
  /**
   * @brief Copies the state of another EnsembleKalmanOutput into this object.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void make_copy(const EnsembleKalmanInversion &another);
  
  // Stored here.
  /** @brief The mcmc. */
  MCMC* mcmc;
  
  //bool sequencer_limit_is_fixed;
  
  // stored here
  /** @brief The index. */
  Index* index;
  
  // If this is not 1.0, then check to see if we think the current ensemble is Gaussian, then skip to the end of the sequence if it is.
  /** @brief The significance level. */
  double significance_level;
  
  /** @brief The estimator type. */
  size_t estimator_type;
  
  //void smc_step();
  
  //void weight_update();
};
}

#endif
