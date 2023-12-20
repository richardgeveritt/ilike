#ifndef KALMANFILTER_H
#define KALMANFILTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"

class KalmanFilterOutput;
class KalmanUpdater;
class KalmanPredictor;

class KalmanFilter : public LikelihoodEstimator
{

public:

  KalmanFilter();

  KalmanFilter(Data* data_in,
               size_t lag_in,
               const std::string &state_name,
               const arma::colvec &prior_mean_in,
               const arma::mat &prior_covariance_in,
               const std::string &index_name_in,
               const std::string &time_name_in,
               const std::string &time_diff_name_in,
               const std::vector<std::string> &measurements_names_in,
               size_t first_index_in,
               size_t last_index_in,
               size_t predictions_per_update_in,
               double update_time_step_in,
               double initial_time_in,
               bool last_index_is_fixed_in,
               KalmanPredictor* predictor_in,
               KalmanUpdater* updater_in,
               bool smcfixed_flag_in,
               const std::string &results_name_in);

  virtual ~KalmanFilter();

  KalmanFilter(const KalmanFilter &another);

  void operator=(const KalmanFilter &another);
  //LikelihoodEstimator* duplicate() const;

  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  LikelihoodEstimator* duplicate() const;
  
  void setup();
  void setup(const Parameters &parameters);
  
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
  
  arma::colvec prior_mean;
  arma::mat prior_covariance;
  
  std::string state_name;
  size_t state_dimension;
  std::string index_name;
  std::string time_name;
  std::string time_diff_name;
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
  
  Parameters schedule_parameters;
  
  std::string results_name;

  void make_copy(const KalmanFilter &another);

  std::ofstream log_likelihood_file_stream;
  std::ofstream time_file_stream;
  std::ofstream vector_variables_file_stream;
  std::ofstream vector_variable_sizes_file_stream;
  std::ofstream output_lengths_file_stream;
  std::ofstream incremental_log_likelihood_file_stream;
  std::ofstream schedule_parameters_file_stream;
  std::ofstream posterior_means_file_stream;
  std::ofstream posterior_covariances_file_stream;
  std::ofstream predicted_means_file_stream;
  std::ofstream predicted_covariances_file_stream;
  
};

#endif
