#ifndef GENERICMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define GENERICMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

class GenericMeasurementCovarianceEstimator;

class GenericMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{

public:

  GenericMeasurementCovarianceEstimatorOutput();
  virtual ~GenericMeasurementCovarianceEstimatorOutput();
  
  GenericMeasurementCovarianceEstimatorOutput(GenericMeasurementCovarianceEstimator* generic_estimator_in);

  GenericMeasurementCovarianceEstimatorOutput(const GenericMeasurementCovarianceEstimatorOutput &another);

  void operator=(const GenericMeasurementCovarianceEstimatorOutput &another);
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  
  void specific_simulate(const Parameters &parameters);
  void subsample_specific_simulate(const Parameters &parameters);
  
  arma::rowvec get_measurement_state_for_covariance() const;
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  arma::colvec get_deterministic_shift() const;
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,const arma::mat &inv_sigma_precomp,
                                            double log_det_precomp) const;
  /*
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const Parameters &conditioned_on_parameters);
  */
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const arma::mat &inv_sigma_precomp,
                                                      double log_det_precomp) const;
  /*
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const Parameters &conditioned_on_parameters);
  */
  
  double evaluate_likelihood();
  /*
  double evaluate_likelihood(const Parameters &conditioned_on_parameters);
  */
  double subsample_evaluate_likelihood();
  /*
  double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);
  */
  
  MeasurementCovarianceEstimator* get_estimator();
  
  void close_ofstreams();

protected:
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  GenericMeasurementCovarianceEstimator* generic_estimator;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);
  
  //Parameters simulated_measurement;
  arma::colvec measurement_state;
  arma::colvec random_shift;

  void make_copy(const GenericMeasurementCovarianceEstimatorOutput &another);
  
};

#endif
