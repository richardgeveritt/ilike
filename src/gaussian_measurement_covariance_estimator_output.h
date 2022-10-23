#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

class GaussianMeasurementCovarianceEstimator;

class GaussianMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{

public:

  GaussianMeasurementCovarianceEstimatorOutput();
  virtual ~GaussianMeasurementCovarianceEstimatorOutput();
  
  GaussianMeasurementCovarianceEstimatorOutput(GaussianMeasurementCovarianceEstimator* gaussian_estimator);

  GaussianMeasurementCovarianceEstimatorOutput(const GaussianMeasurementCovarianceEstimatorOutput &another);

  void operator=(const GaussianMeasurementCovarianceEstimatorOutput &another);
  virtual GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const=0;
  
  arma::rowvec get_measurement_state_for_covariance() const;
  
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  
  arma::colvec get_deterministic_shift() const;
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Ginv,
                           const arma::mat &Ftranspose,
                           const arma::mat &V,
                           double inverse_incremental_temperature);
  
  arma::mat get_kalman_gain(const arma::mat &Cxy,
                            const arma::mat &Cyy,
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
  
  virtual arma::mat get_measurement_covariance()=0;

protected:
  
  arma::colvec measurement_state;
  arma::colvec random_shift;
  //arma::mat measurement_noise;

  void make_copy(const GaussianMeasurementCovarianceEstimatorOutput &another);
  
};

#endif
