#ifndef ENSEMBLEKALMANMFDS_H
#define ENSEMBLEKALMANMFDS_H

#include "ensemble_kalman.h"
#include "ensemble_sequencer.h"

namespace ilike
{
  /**
   * @file ensemble_kalman_mfds.h
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
class Index;

class EnsembleKalmanMFDS : public EnsembleKalman
{
public:
  
  /**
   * @brief Performs the ensemblekalmanmfds operation.
   */
  EnsembleKalmanMFDS();
  
  // constructor for all variants that have a joint representation as intractable and Gaussian
  EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                     size_t* seed_in,
                     Data* data_in,
                     size_t lag_in,
                     const std::vector<Parameters> &initial_points_in,
                     EnsembleShifter* shifter_in,
                     double delta_t_in,
                     size_t number_of_iterations_in,
                     const Parameters &prior_means,
                     const Parameters &prior_covariances,
                     SimulateModelPtr simulate_model_in,
                     size_t update_type,
                     std::shared_ptr<Transform> summary_statistics_in,
                     std::shared_ptr<Transform> transform_in,
                     const std::vector<std::string> &measurement_variables_in,
                     bool parallel_in,
                     size_t grain_size_in,
                     const std::string &results_name_in);
  
  // constructor for case of everything Gaussians - must construct data and Gaussian info in the correct way (data must include prior, etc)
  EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                     size_t* seed_in,
                     Data* data_in,
                     size_t lag_in,
                     const std::vector<Parameters> &initial_points_in,
                     EnsembleShifter* shifter_in,
                     double delta_t_in,
                     size_t number_of_iterations_in,
                     std::shared_ptr<Transform> measurement_transform_function_in,
                     const std::vector<std::string> &measurement_variables,
                     const std::vector<GetMatrixPtr> &measurement_noise_functions,
                     std::shared_ptr<Transform> summary_statistics_in,
                     std::shared_ptr<Transform> transform_in,
                     bool parallel_in,
                     size_t grain_size_in,
                     const std::string &results_name_in);
  
  // constructor for case of everything intractable - must construct data and projection to measurements in the correct way (data must include prior, etc)
  EnsembleKalmanMFDS(RandomNumberGenerator* rng_in,
                     size_t* seed_in,
                     Data* data_in,
                     size_t lag_in,
                     const std::vector<Parameters> &initial_points_in,
                     EnsembleShifter* shifter_in,
                     double delta_t_in,
                     size_t number_of_iterations_in,
                     SimulateModelPtr simulate_model_in,
                     std::shared_ptr<Transform> summary_statistics_in,
                     std::shared_ptr<Transform> transform_in,
                     const std::vector<std::string> &measurement_variables_in,
                     bool parallel_in,
                     size_t grain_size_in,
                     const std::string &results_name_in);
  
  /**
   * @brief Performs the ensemblekalmanmfds operation.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  EnsembleKalmanMFDS(const EnsembleKalmanMFDS &another);
  
  /**
   * @brief Performs the ~ensemblekalmanmfds operation.
   *
   * @param void The void.
   */
  virtual ~EnsembleKalmanMFDS(void);
  
  /**
   * @brief Assignment operator for EnsembleKalmanOutput.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void operator=(const EnsembleKalmanMFDS &another);
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
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle);
  
  /*
   MoveOutput* move(RandomNumberGenerator &rng,
   Particle &particle,
   const Parameters &conditioned_on_parameters);
   */
  
  /*
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
  
  /**
   * @brief Performs the ensemble kalman subsample simulate operation.
   *
   * @param simulation The simulation.
   */
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                          const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the predict operation.
   *
   * @param simulation The simulation.
   */
  void predict(EnsembleKalmanOutput* simulation);
  
  //void setup_variables();
  
  //void smc_update(EnsembleKalmanOutput* current_state);
  
  /**
   * @brief Copies the state of another EnsembleKalmanOutput into this object.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void make_copy(const EnsembleKalmanMFDS &another);
  
  // Stored here.
  /** @brief The mcmc. */
  MCMC* mcmc;
  
  /** @brief The state mean as data. */
  Data state_mean_as_data;
  
  /** @brief The sequencer. */
  EnsembleSequencer sequencer;
  /** @brief The sequencer limit is fixed. */
  bool sequencer_limit_is_fixed;
  
  /** @brief The delta t. */
  double delta_t;
  /** @brief The number of iterations. */
  size_t number_of_iterations;
  
  //arma::colvec predicted_mean;
  
  // stored here
  /** @brief The index. */
  Index* index;
  
  //void smc_step();
  
  //void weight_update();
};
}

#endif
