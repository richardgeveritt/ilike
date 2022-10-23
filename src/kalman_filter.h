#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "likelihood_estimator.h"
#include "function_pointers.h"
#include "parameters.h"

class KalmanFilterOutput;
class KalmanUpdater;
class KalmanPredictor;

class KalmanFilter : public LikelihoodEstimator
{

public:

  KalmanFilter();

  KalmanFilter(RandomNumberGenerator* rng_in,
               size_t* seed_in,
               Data* data_in,
               size_t lag_in,
               EvaluateLogLikelihoodPtr llhd_in,
               double current_time_in=0.0,
               bool smcfixed_flag_in = TRUE);

  virtual ~KalmanFilter();

  KalmanFilter(const KalmanFilter &another);

  void operator=(const KalmanFilter &another);
  //LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  LikelihoodEstimatorOutput* initialise();
  KalmanFilterOutput* kalman_filter_initialise();
  
  void evaluate(KalmanFilterOutput* simulation);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  KalmanFilterOutput* run();

  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  KalmanFilterOutput* kalman_filter_initialise(const Parameters &parameters);
  
  void evaluate(KalmanFilterOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  void subsample_evaluate(KalmanFilterOutput* simulation);
  
  void subsample_evaluate(KalmanFilterOutput* simulation,
                const Parameters &conditioned_on_parameters);

  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  KalmanFilterOutput* run(const Parameters &conditioned_on_parameters);

protected:
  
  bool check_termination() const;

  friend KalmanFilterOutput;
  
  // A flag to determine if the terms are deemed to be "smcfixed" (will not be reevaluated when finding the next target in adaptive SMC).
  bool smcfixed_flag;
  
  arma::colvec prior_mean;
  arma::mat prior_covariance;
  
  std::string index_name;
  //std::string time_name;
  std::vector<std::string> measurements_names;
  size_t first_index;
  size_t last_index;
  size_t predictions_per_update;
  double update_time_step;
  double current_time;
  size_t current_index;
  bool last_index_is_fixed;
  size_t lag;
  
  // stored here
  KalmanUpdater* updater;
  KalmanPredictor* predictor;

  void make_copy(const KalmanFilter &another);

};

#endif
