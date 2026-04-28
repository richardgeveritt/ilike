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
  /**
   * @file kalman_filter.h
   * @brief Defines the KalmanFilterOutput class.
   *
   * Stores and manages the output produced by KalmanFilter. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class KalmanFilterOutput
   * @brief The kalman filter output class.
   */


class KalmanFilterOutput;
class KalmanUpdater;
class KalmanPredictor;

class KalmanFilter : public LikelihoodEstimator
{
  
public:
  
  /**
   * @brief Performs the kalmanfilter operation.
   */
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
  
  /**
   * @brief Performs the ~kalmanfilter operation.
   */
  virtual ~KalmanFilter();
  
  /**
   * @brief Performs the kalmanfilter operation.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  KalmanFilter(const KalmanFilter &another);
  
  /**
   * @brief Assignment operator for KalmanFilterOutput.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void operator=(const KalmanFilter &another);
  //LikelihoodEstimator* duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  /**
   * @brief Creates a deep copy of this KalmanFilterOutput object.
   *
   * @return The result.
   */
  LikelihoodEstimator* duplicate() const;
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  void setup();
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  void setup(const Parameters &parameters);
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise();
  /**
   * @brief Performs the kalman filter initialise operation.
   *
   * @return The result.
   */
  KalmanFilterOutput* kalman_filter_initialise();
  
  /**
   * @brief Evaluates.
   *
   * @param simulation The simulation.
   */
  void evaluate(KalmanFilterOutput* simulation);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @return The result.
   */
  KalmanFilterOutput* run();
  
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  /**
   * @brief Performs the kalman filter initialise operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  KalmanFilterOutput* kalman_filter_initialise(const Parameters &parameters);
  
  void evaluate(KalmanFilterOutput* simulation,
                const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the subsample evaluate operation.
   *
   * @param simulation The simulation.
   */
  void subsample_evaluate(KalmanFilterOutput* simulation);
  
  void subsample_evaluate(KalmanFilterOutput* simulation,
                          const Parameters &conditioned_on_parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
  /**
   * @brief Runs the algorithm and returns the collected output.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  KalmanFilterOutput* run(const Parameters &conditioned_on_parameters);
  
protected:
  
  /**
   * @brief Class-specific implementation for change data.
   *
   * @param new_data The new data.
   */
  void specific_change_data(Data* new_data);
  
  /**
   * @brief Performs the check termination operation.
   *
   * @return The result.
   */
  bool check_termination() const;
  
  friend KalmanFilterOutput;
  
  /** @brief The prior mean. */
  arma::colvec prior_mean;
  /** @brief The prior covariance. */
  arma::mat prior_covariance;
  
  /** @brief The state name. */
  std::string state_name;
  /** @brief The state dimension. */
  size_t state_dimension;
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
  /** @brief The last index is fixed. */
  bool last_index_is_fixed;
  /** @brief The lag. */
  size_t lag;
  
  // stored here
  /** @brief The updaters. */
  std::vector<KalmanUpdater*> updaters;
  /** @brief The predictor. */
  KalmanPredictor* predictor;
  
  /** @brief The schedule parameters. */
  Parameters schedule_parameters;
  
  /** @brief The results name. */
  std::string results_name;
  
  /**
   * @brief Copies the state of another KalmanFilterOutput into this object.
   *
   * @param another The KalmanFilterOutput instance to copy from.
   */
  void make_copy(const KalmanFilter &another);
  
  /** @brief HDF5 output file (kept open for the duration of a run). */
  std::shared_ptr<HighFive::File> h5_file;
  /** @brief Path to the HDF5 output file. */
  std::string h5_file_path;
  
};
}

#endif
