#ifndef MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"

class Parameters;
class MeasurementCovarianceEstimator;
class Index;

class MeasurementCovarianceEstimatorOutput
{

public:

  MeasurementCovarianceEstimatorOutput();
  MeasurementCovarianceEstimatorOutput(MeasurementCovarianceEstimator* estimator_in);
  virtual ~MeasurementCovarianceEstimatorOutput();

  MeasurementCovarianceEstimatorOutput(const MeasurementCovarianceEstimatorOutput &another);

  void operator=(const MeasurementCovarianceEstimatorOutput &another);
  virtual MeasurementCovarianceEstimatorOutput* duplicate() const=0;
  
  virtual void simulate(const Parameters &parameters)=0;
  virtual arma::rowvec get_measurement_state_for_covariance() const=0;
  virtual arma::colvec get_shift(double inverse_incremental_temperature) const=0;
  virtual arma::colvec get_deterministic_shift() const=0;
  virtual arma::mat get_kalman_gain(const arma::mat &Cxy,
                                    const arma::mat &Cyy,
                                    double inverse_incremental_temperature)=0;
  virtual arma::mat get_adjustment(const arma::mat &Zf,
                                   const arma::mat &Ginv,
                                   const arma::mat &Ftranspose,
                                   const arma::mat &V,
                                   double inverse_incremental_temperature)=0;
  arma::colvec* get_measurement() const;
  
  virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature)=0;
  virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                    const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                              const Parameters &conditioned_on_parameters)=0;
  
  virtual double evaluate_likelihood()=0;
  virtual double evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_likelihood()=0;
  virtual double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
protected:
  
  // Stored in ModelAndAlgorithm or in main.
  MeasurementCovarianceEstimator* estimator;
  
  void make_copy(const MeasurementCovarianceEstimatorOutput &another);
  
  bool set_using_parameters;

};

#endif
