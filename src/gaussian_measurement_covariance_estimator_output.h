#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator_output.h"

namespace ilike
{
class GaussianMeasurementCovarianceEstimator;

class GaussianMeasurementCovarianceEstimatorOutput : public MeasurementCovarianceEstimatorOutput
{
  
public:
  
  GaussianMeasurementCovarianceEstimatorOutput();
  virtual ~GaussianMeasurementCovarianceEstimatorOutput();
  
  GaussianMeasurementCovarianceEstimatorOutput(const GaussianMeasurementCovarianceEstimatorOutput &another);
  
  void operator=(const GaussianMeasurementCovarianceEstimatorOutput &another);
  virtual GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const=0;
  
  arma::rowvec get_measurement_state_for_covariance() const;
  
  arma::colvec get_shift(double inverse_incremental_temperature) const;
  
  arma::colvec get_deterministic_shift() const;
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature);
  
  double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                            const arma::mat &inv_sigma_precomp,
                                            double log_det_precomp) const;
  
  double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                      const arma::mat &inv_sigma_precomp,
                                                      double log_det_precomp) const;
  
  /*
   double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters);
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
  
  //virtual arma::mat get_measurement_covariance()=0;
  
  virtual GaussianMeasurementCovarianceEstimator* get_gaussian_estimator()=0;
  virtual const GaussianMeasurementCovarianceEstimator* get_gaussian_estimator() const=0;
  
  void close_ofstreams();
  
protected:
  
  arma::colvec measurement_state;
  arma::colvec random_shift;
  
  //arma::mat measurement_noise;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index);
  
  void make_copy(const GaussianMeasurementCovarianceEstimatorOutput &another);
  
};
}

#endif
