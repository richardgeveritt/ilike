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
  
  void simulate(const Parameters &parameters);
  
  arma::rowvec get_measurement_state_for_covariance() const;
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  arma::colvec get_deterministic_shift() const;
  
  arma::mat get_kalman_gain(const arma::mat &Cxy,
                            const arma::mat &Cyy,
                            double inverse_incremental_temperature);
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Ginv,
                           const arma::mat &Ftranspose,
                           const arma::mat &V,
                           double inverse_incremental_temperature);
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature);
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const Parameters &conditioned_on_parameters);
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const Parameters &conditioned_on_parameters);
  
  double evaluate_likelihood();
  double evaluate_likelihood(const Parameters &conditioned_on_parameters);
  double subsample_evaluate_likelihood();
  double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters);

protected:
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  GenericMeasurementCovarianceEstimator* generic_estimator;
  
  //Parameters simulated_measurement;
  arma::colvec measurement_state;
  arma::colvec random_shift;

  void make_copy(const GenericMeasurementCovarianceEstimatorOutput &another);
  
};

#endif
