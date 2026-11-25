#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include <memory>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"
#include "ilike_hdf5_utils.h"

namespace ilike
{

class KalmanFilterOutput;
class KalmanUpdater;
class KalmanPredictor;

class KalmanFilter : public LikelihoodEstimator
{
  
public:
  
  KalmanFilter();
  
  KalmanFilter(Data* data_in,
               size_t lag_in,
               //const std::string &state_name,
               const arma::colvec &prior_mean_in,
               const arma::mat &prior_covariance_in,
               const std::string &index_name_in,
               const std::string &time_name_in,
               const std::string &time_diff_name_in,
               //const std::vector<std::string> &measurements_names_in,
               size_t first_index_in,
               size_t last_index_in,
               size_t predictions_per_update_in,
               double update_time_step_in,
               double initial_time_in,
               bool last_index_is_fixed_in,
               KalmanPredictor* predictor_in,
               const std::vector<KalmanUpdater*> &updaters_in,
               bool smcfixed_flag_in,
               const std::string &results_name_in);
  
  virtual ~KalmanFilter();
  
  KalmanFilter(const KalmanFilter &another);
  
  void operator=(const KalmanFilter &another);

  LikelihoodEstimator* duplicate() const;
  
  void setup();
  void setup(const Parameters &parameters);
  
  LikelihoodEstimatorOutput* initialise();
  KalmanFilterOutput* kalman_filter_initialise();
  
  void evaluate(KalmanFilterOutput* simulation);

  KalmanFilterOutput* run();
  
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  KalmanFilterOutput* kalman_filter_initialise(const Parameters &parameters);
  
  void evaluate(KalmanFilterOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(KalmanFilterOutput* simulation);
  
  void subsample_evaluate(KalmanFilterOutput* simulation,
                          const Parameters &conditioned_on_parameters);

  KalmanFilterOutput* run(const Parameters &conditioned_on_parameters);

  // (dev/plan.md 3.21 (c2''')) SMC^2 with a Kalman filter: the caller's data
  // tempering takes in the filter's observations one at a time (or a share of
  // one, with anneal). Each output keeps its own committed statistics; a trial
  // at the caller's target runs one predict-and-update step on a copy of the
  // output (no random numbers: the increment is the exact predictive density);
  // the caller's move stage commits the observation once it is whole; a new
  // output is run through the observations committed so far.
  void set_follows_caller_target(const std::string &fraction_variable, bool anneal);
  bool advances_at_caller_target() const { return this->follows_caller_target; }
  // (C270) parallel outer chunks: a copy (its own predictor and updaters,
  // and the shared state an evaluation changes) runs in each chunk group; the
  // filter draws no random numbers; its outputs are repointed afterwards
  bool supports_parallel_copies() const { return true; }
  size_t observations() const { return this->last_index - this->first_index + 1; }
  void advance_at_caller_target_kf(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters);
  
protected:

  bool follows_caller_target = false;
  std::string followed_variable;
  bool followed_anneal = false;
  size_t followed_steps = 0; // (the most observations any output has committed)
  void follow_evaluate(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters);
  void follow_step(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters);
  double follow_trial_value(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters,
                            size_t k, double r);
  void follow_commit(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters);
  void follow_from_start(KalmanFilterOutput* current_state, const Parameters &conditioned_on_parameters);
  
  void specific_change_data(Data* new_data);
  
  bool check_termination() const;
  
  friend KalmanFilterOutput;
  
  arma::colvec prior_mean;
  arma::mat prior_covariance;
  
  std::string state_name;
  size_t state_dimension;
  std::string index_name;
  std::string time_name;
  std::string time_diff_name;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;
  bool last_index_is_fixed;
  size_t lag;
  
  // stored here
  // owned (dev/plan.md 3.5, CH3)
  std::vector<std::unique_ptr<KalmanUpdater>> updaters;
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<KalmanPredictor> predictor;
  
  Parameters schedule_parameters;
  
  std::string results_name;
  
  void make_copy(const KalmanFilter &another);
  
  std::shared_ptr<HighFive::File> h5_file;
  std::string h5_file_path;
  
};
}

#endif
