#ifndef ITERATIVEENSEMBLEKALMANINVERSION_H
#define ITERATIVEENSEMBLEKALMANINVERSION_H

#include "ensemble_kalman.h"
#include "ensemble_sequencer.h"

class EnsembleKalmanOutput;
class MCMC;
class MoveOutput;

class EnsembleKalmanInversion : public EnsembleKalman
{
public:

  EnsembleKalmanInversion();
  
  // single move from start to end
  EnsembleKalmanInversion(RandomNumberGenerator* rng_in,
                          size_t* seed_in,
                          Data* data_in,
                          size_t number_of_ensemble_members_in,
                          size_t lag_in,
                          EnsembleShifter* shifter_in,
                          SimulateDistributionPtr simulate_prior_in,
                          SimulateModelPtr simulate_model_in,
                          std::shared_ptr<Transform> transform_in,
                          std::shared_ptr<Transform> summary_statistics_in,
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
                          SimulateDistributionPtr simulate_prior_in,
                          SimulateModelPtr simulator_model_in,
                          std::shared_ptr<Transform> transform_in,
                          std::shared_ptr<Transform> summary_statistics_in,
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
                          SimulateDistributionPtr simulate_prior_in,
                          std::shared_ptr<Transform> measurement_transform_function_in,
                          const std::vector<std::string> &measurement_variables,
                          const std::vector<GetMeasurementMatrixPtr> &measurement_noise_functions,
                          std::shared_ptr<Transform> transform_in,
                          std::shared_ptr<Transform> summary_statistics_in,
                          bool parallel_in,
                          size_t grain_size_in,
                          const std::string &results_name_in);
  
  EnsembleKalmanInversion(const EnsembleKalmanInversion &another);
  
  virtual ~EnsembleKalmanInversion(void);

  void operator=(const EnsembleKalmanInversion &another);
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
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
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

  void make_copy(const EnsembleKalmanInversion &another);
  
  // Stored here.
  MCMC* mcmc;
  
  //bool sequencer_limit_is_fixed;
  
  // stored here
  Index* index;

  //void smc_step();

  //void weight_update();
};

#endif
