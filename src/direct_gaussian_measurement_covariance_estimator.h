#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"
#include "function_pointers.h"
#include "gaussian_independent_proposal_kernel.h"

class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectGaussianMeasurementCovarianceEstimator : public GaussianMeasurementCovarianceEstimator
{

public:

  DirectGaussianMeasurementCovarianceEstimator();

  virtual ~DirectGaussianMeasurementCovarianceEstimator();

  DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another);

  void operator=(const DirectGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;

protected:
  
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);

  void make_copy(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  TransformPtr transform_function;
  std::vector<GetMeasurementMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;

};

#endif
