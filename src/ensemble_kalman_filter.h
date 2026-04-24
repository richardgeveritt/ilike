#ifndef ENSEMBLEKALMANFILTER_H
#define ENSEMBLEKALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "ensemble_kalman.h"
#include "ilike_header.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file ensemble_kalman_filter.h
   * @brief Defines the EnsembleKalmanOutput class.
   *
   * Stores and manages the output produced by EnsembleKalman. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class EnsembleKalmanOutput
   * @brief The ensemble kalman output class.
   */


class EnsembleKalmanOutput;
class ProposalKernel;
class MoveOutput;
//class EnsembleKalmanUpdater;
//class EnsembleKalmanPredictor;

class EnsembleKalmanFilter : public EnsembleKalman
{
  
public:
  
  /**
   * @brief Performs the ensemblekalmanfilter operation.
   */
  EnsembleKalmanFilter();
  
  EnsembleKalmanFilter(RandomNumberGenerator* rng_in,
                       size_t* seed_in,
                       Data* data_in,
                       size_t lag_in,
                       const std::string &index_name_in,
                       const std::string &time_name_in,
                       const std::string &time_diff_name_in,
                       size_t first_index_in,
                       size_t last_index_in,
                       size_t predictions_per_update_in,
                       double update_time_step_in,
                       double initial_time_in,
                       size_t number_of_ensemble_members_in,
                       EnsembleShifter* shifter_in,
                       std::shared_ptr<Transform> transform_in,
                       IndependentProposalKernel* prior_in,
                       ProposalKernel* transition_model_in,
                       const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in,
                       bool smcfixed_flag_in,
                       bool sequencer_limit_is_fixed_in,
                       bool parallel_in,
                       size_t grain_size_in,
                       const std::string &results_name_in);
  
  /**
   * @brief Performs the ~ensemblekalmanfilter operation.
   */
  virtual ~EnsembleKalmanFilter();
  
  /**
   * @brief Performs the ensemblekalmanfilter operation.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  EnsembleKalmanFilter(const EnsembleKalmanFilter &another);
  
  /**
   * @brief Assignment operator for EnsembleKalmanOutput.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void operator=(const EnsembleKalmanFilter &another);
  /**
   * @brief Creates a deep copy of this EnsembleKalmanOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a ensemble_kalman pointer.
   *
   * @return The result.
   */
  EnsembleKalman* ensemble_kalman_duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise();
  /**
   * @brief Class-specific implementation for ensemble kalman initialise.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters);
  
  //void evaluate(EnsembleKalmanOutput* simulation);
  
  void evaluate(EnsembleKalmanOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(EnsembleKalmanOutput* simulation,
                          const Parameters &conditioned_on_parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  /**
   * @brief Class-specific implementation for run.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_run();
  /**
   * @brief Class-specific implementation for run.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters);
  
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
  
  //void weight_for_adapting_sequence(Ensemble &current_particles,
  //                                  double incremental_temperature);
  
protected:
  
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
  
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
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
  //void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                          const Parameters &conditioned_on_parameters);
  
  void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  /*
   void ensemble_kalman_subsample_simulate(EnsembleKalmanOutput* simulation,
   const Parameters &conditioned_on_parameters);
   void ensemble_kalman_subsample_evaluate(EnsembleKalmanOutput* simulation,
   const Parameters &conditioned_on_parameters);
   
   void ensemble_kalman_subsample_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
   const Parameters &conditioned_on_parameters);
   void ensemble_kalman_subsample_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
   const Parameters &conditioned_on_parameters);
   */
  
  bool check_termination() const;
  
  //
  /**
   * @brief Performs the setup variables operation.
   */
  void setup_variables();
  
  friend EnsembleKalmanOutput;
  
  
  //Parameters schedule_parameters;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  /** @brief The smcfixed flag. */
  bool smcfixed_flag;
  
  // point to the one in the hmm_ensemble_factors required to construct this class
  /** @brief The transition model kernel. */
  ProposalKernel* transition_model_kernel;
  
  // stored here
  /** @brief The index. */
  Index* index;
  
  //arma::colvec prior_mean;
  //arma::mat prior_covariance;
  
  /** @brief The index name. */
  std::string index_name;
  /** @brief The time name. */
  std::string time_name;
  /** @brief The time diff name. */
  std::string time_diff_name;
  //std::vector<std::string> measurements_names;
  /** @brief The first index. */
  size_t first_index;
  /** @brief The last index. */
  size_t last_index;
  /** @brief The predictions per update. */
  size_t predictions_per_update;
  /** @brief The update time step. */
  double update_time_step;
  /** @brief The current time. */
  double current_time;
  /** @brief The current index. */
  size_t current_index;
  //bool last_index_is_fixed;
  
  
  // Stored here.
  //EnsembleKalmanFilterOutput* output;
  
  /**
   * @brief Copies the state of another EnsembleKalmanOutput into this object.
   *
   * @param another The EnsembleKalmanOutput instance to copy from.
   */
  void make_copy(const EnsembleKalmanFilter &another);
  
};
}

#endif
