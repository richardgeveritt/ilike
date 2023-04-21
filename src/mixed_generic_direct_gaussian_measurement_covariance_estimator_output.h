#ifndef MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

class MixedGenericDirectGaussianMeasurementCovarianceEstimator;

class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{

public:

  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput();
  virtual ~MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput();
  
  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(MixedGenericDirectGaussianMeasurementCovarianceEstimator* estimator_in);

  MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);

  void operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  
  void specific_simulate(const Parameters &parameters);
  void subsample_specific_simulate(const Parameters &parameters);
  
  arma::rowvec get_measurement_state_for_covariance() const;
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  arma::colvec get_deterministic_shift() const;
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature);
  /*
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const Parameters &conditioned_on_parameters);
  */
  
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature);
  
  /*
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const Parameters &conditioned_on_parameters);
  */
  
  double evaluate_likelihood();
  //double evaluate_likelihood(const Parameters &conditioned_on_parameters);
  double subsample_evaluate_likelihood();
  //double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);
  
  MeasurementCovarianceEstimator* get_estimator();
  
  void close_ofstreams();

protected:
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  MixedGenericDirectGaussianMeasurementCovarianceEstimator* estimator;
  
  //Parameters simulated_measurement;
  arma::colvec likelihood_measurement_state;
  arma::colvec prior_measurement_state;
  arma::colvec likelihood_random_shift;
  arma::colvec prior_random_shift;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);

  void make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput &another);
  
};

#endif
