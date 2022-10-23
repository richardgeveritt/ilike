#ifndef ENSEMBLEKALMANFILTER_H
#define ENSEMBLEKALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "ensemble_kalman.h"
#include "function_pointers.h"
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
                       EvaluateLogLikelihoodPtr llhd_in,
                       double current_time_in,
                       bool sequencer_limit_is_fixed_in);

  virtual ~EnsembleKalmanFilter();

  EnsembleKalmanFilter(const EnsembleKalmanFilter &another);

  void operator=(const EnsembleKalmanFilter &another);
  //LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;

  EnsembleKalmanOutput* ensemble_kalman_initialise();
  EnsembleKalmanOutput* ensemble_kalman_initialise(const Parameters &parameters);
  
  /*
  void evaluate(EnsembleKalmanOutput* simulation);
  void evaluate(EnsembleKalmanOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(EnsembleKalmanOutput* simulation,
                          const Parameters &conditioned_on_parameters);
  */

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  EnsembleKalmanOutput* specific_run();
  EnsembleKalmanOutput* specific_run(const Parameters &conditioned_on_parameters);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle);
  
  MoveOutput* move(RandomNumberGenerator &rng,
                   Particle &particle,
                   const Parameters &conditioned_on_parameters);
  
  MoveOutput* subsample_move(RandomNumberGenerator &rng,
                             Particle &particle,
                             const Parameters &conditioned_on_parameters);
  
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
  
  void ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                              const Parameters &conditioned_on_parameters);
  void ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                const Parameters &conditioned_on_parameters);
  
  void subsample_ensemble_kalman_simulate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  void subsample_ensemble_kalman_evaluate(EnsembleKalmanOutput* simulation,
                                          const Parameters &conditioned_on_parameters);
  
  void subsample_ensemble_kalman_evaluate_smcfixed_part(EnsembleKalmanOutput* simulation,
                                                        const Parameters &conditioned_on_parameters);
  void subsample_ensemble_kalman_evaluate_smcadaptive_part_given_smcfixed(EnsembleKalmanOutput* simulation,
                                                                          const Parameters &conditioned_on_parameters);
  
  bool check_termination() const;

  friend EnsembleKalmanOutput;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  // point to the one in the hmm_ensemble_factors required to construct this class
  ProposalKernel* proposal_kernel;
  
  //arma::colvec prior_mean;
  //arma::mat prior_covariance;
  
  std::string index_name;
  //std::string time_name;
  //std::vector<std::string> measurements_names;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;
  bool last_index_is_fixed;


  // Stored here.
  //EnsembleKalmanFilterOutput* output;

  void make_copy(const EnsembleKalmanFilter &another);

};

#endif
