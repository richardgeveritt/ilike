#ifndef ENSEMBLEKALMANMFDS_H
#define ENSEMBLEKALMANMFDS_H

#include "ensemble_kalman.h"
#include "ensemble_sequencer.h"

namespace ilike
{
class EnsembleKalmanOutput;
class MCMC;
class Index;

class EnsembleKalmanMFDS : public EnsembleKalman
{
public:
  
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
  
  EnsembleKalmanMFDS(const EnsembleKalmanMFDS &another);
  
  virtual ~EnsembleKalmanMFDS(void);
  
  void operator=(const EnsembleKalmanMFDS &another);
  EnsembleKalman* ensemble_kalman_duplicate() const;
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
  
  EnsembleKalmanOutput* specific_run();
  
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise();
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation);
  
  EnsembleKalmanOutput* specific_run(const Parameters &parameters);
  
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters);
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
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
  
  void predict(EnsembleKalmanOutput* simulation);
  
  //void setup_variables();
  
  //void smc_update(EnsembleKalmanOutput* current_state);
  
  void make_copy(const EnsembleKalmanMFDS &another);
  
  // Stored here.
  MCMC* mcmc;
  
  Data state_mean_as_data;
  
  EnsembleSequencer sequencer;
  bool sequencer_limit_is_fixed;
  
  double delta_t;
  size_t number_of_iterations;
  
  //arma::colvec predicted_mean;
  
  // stored here
  Index* index;
  
  //void smc_step();
  
  //void weight_update();
};
}

#endif
