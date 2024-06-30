#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

//#include <RcppArmadillo.h>
//using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator_output.h"
#include "gaussian_independent_proposal_kernel.h"

class DirectGaussianMeasurementCovarianceEstimator;

class DirectGaussianMeasurementCovarianceEstimatorOutput : public GaussianMeasurementCovarianceEstimatorOutput
{

public:

  DirectGaussianMeasurementCovarianceEstimatorOutput();
  
  DirectGaussianMeasurementCovarianceEstimatorOutput(DirectGaussianMeasurementCovarianceEstimator* direct_estimator_in);

  virtual ~DirectGaussianMeasurementCovarianceEstimatorOutput();

  DirectGaussianMeasurementCovarianceEstimatorOutput(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);

  void operator=(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);
  MeasurementCovarianceEstimatorOutput* duplicate() const;
  GaussianMeasurementCovarianceEstimatorOutput* gaussian_duplicate() const;
  
  void specific_simulate(const Parameters &parameters);
  void subsample_specific_simulate(const Parameters &parameters);
  
  //arma::rowvec get_measurement_random_shift();
  
  //arma::mat get_measurement_covariance() const;
  
  MeasurementCovarianceEstimator* get_estimator();
  GaussianMeasurementCovarianceEstimator* get_gaussian_estimator();
  const GaussianMeasurementCovarianceEstimator* get_gaussian_estimator() const;

protected:
  
  // not stored here
  DirectGaussianMeasurementCovarianceEstimator* direct_estimator;

  void make_copy(const DirectGaussianMeasurementCovarianceEstimatorOutput &another);

};

#endif
