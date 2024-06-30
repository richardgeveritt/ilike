#ifndef ENSEMBLEKALMANFILTER_H
#define ENSEMBLEKALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "ensemble_kalman.h"
#include "ilike_header.h"
#include "parameters.h"

class EnsembleKalmanOutput;
class ProposalKernel;
class MoveOutput;
//class EnsembleKalmanUpdater;
//class EnsembleKalmanPredictor;

class EnsembleKalmanFilter : public EnsembleKalman
{

public:

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

  virtual ~EnsembleKalmanFilter();

  EnsembleKalmanFilter(const EnsembleKalmanFilter &another);

  void operator=(const EnsembleKalmanFilter &another);
  LikelihoodEstimator* duplicate() const;
  EnsembleKalman* ensemble_kalman_duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  EnsembleKalmanOutput* specific_ensemble_kalman_initialise();
  EnsembleKalmanOutput* specific_ensemble_kalman_initialise(const Parameters &parameters);
  
  //void evaluate(EnsembleKalmanOutput* simulation);
  
  void evaluate(EnsembleKalmanOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(EnsembleKalmanOutput* simulation,
                          const Parameters &conditioned_on_parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  EnsembleKalmanOutput* specific_run();
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
  
  void ensemble_kalman_simulate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation);
  
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation);
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
  void setup_variables();

  friend EnsembleKalmanOutput;
  
  
  //Parameters schedule_parameters;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  // point to the one in the hmm_ensemble_factors required to construct this class
  ProposalKernel* transition_model_kernel;
  
  // stored here
  Index* index;
  
  //arma::colvec prior_mean;
  //arma::mat prior_covariance;
  
  std::string index_name;
  std::string time_name;
  std::string time_diff_name;
  //std::vector<std::string> measurements_names;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;
  //bool last_index_is_fixed;


  // Stored here.
  //EnsembleKalmanFilterOutput* output;

  void make_copy(const EnsembleKalmanFilter &another);

};

#endif
