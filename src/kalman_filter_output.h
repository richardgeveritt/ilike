#ifndef KALMANFILTEROUTPUT_H
#define KALMANFILTEROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"

namespace ilike
{
class KalmanFilter;

class KalmanFilterOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  KalmanFilterOutput();
  KalmanFilterOutput(KalmanFilter* estimator_in,
                     size_t lag_in,
                     const std::string &results_name_in);
  virtual ~KalmanFilterOutput();
  
  KalmanFilterOutput(const KalmanFilterOutput &another);
  void operator=(const KalmanFilterOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  
  void simulate();
  
  void simulate(const Parameters &parameters);
  //double evaluate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &conditioned_on_parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &conditioned_on_parameters);
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void set_current_predicted_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  void add_predicted_statistics();
  
  void set_current_posterior_statistics(const arma::colvec &latest_mean,
                                        const arma::mat &latest_covariance);
  void add_posterior_statistics();
  
  void add_schedule_parameters(const Parameters &current_schedule_parameters);
  
  void add_log_normalising_constant_ratio(double log_normalising_constant_ratio);
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  void set_time();
  void set_llhd(double llhd_in);
  
  arma::colvec predicted_mean_back() const;
  arma::colvec posterior_mean_back() const;
  arma::mat predicted_covariance_back() const;
  arma::mat posterior_covariance_back() const;
  
  void set_current_predicted_to_be_current_posterior();
  
  size_t predicted_size() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);
  
  void forget_you_were_already_written_to_file();
  
  void terminate();
  
  void increment_kf_iteration();
  
  void close_ofstreams();
  
  void print(std::ostream &os) const;
  
  std::deque<double> times;
  
  std::deque<double> llhds;
  
  std::string results_name;
  
  std::chrono::high_resolution_clock::time_point start_time;
  
protected:
  
  // Stored in ModelAndAlgorithm.
  KalmanFilter* estimator;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index="");
  
  void make_copy(const KalmanFilterOutput &another);
  
  std::deque<arma::colvec> predicted_means;
  std::deque<arma::mat> predicted_covariances;
  std::deque<arma::colvec> posterior_means;
  std::deque<arma::mat> posterior_covariances;
  
  std::deque<Parameters> schedule_parameters;
  
  std::deque<double> log_normalising_constant_ratios;
  
  arma::colvec current_predicted_mean;
  arma::mat current_predicted_covariance;
  arma::colvec current_posterior_mean;
  arma::mat current_posterior_covariance;
  
  size_t lag;
  size_t output_lag;
  
  size_t kf_iteration;
  
  int iteration_written_to_file;
  
};
}

#endif
