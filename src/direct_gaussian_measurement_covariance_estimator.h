#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectGaussianMeasurementCovarianceEstimator : public GaussianMeasurementCovarianceEstimator
{

public:

  DirectGaussianMeasurementCovarianceEstimator();
  
  DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               std::shared_ptr<Transform> transform_in,
                                               std::shared_ptr<Transform> summary_statistics_in,
                                               std::shared_ptr<Transform> transform_function_in,
                                               const std::vector<std::string> &measurement_variables_in,
                                               const std::vector<arma::mat> &measurement_noise_in);
  
  DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                               size_t* seed_in,
                                               Data* data_in,
                                               std::shared_ptr<Transform> transform_in,
                                               std::shared_ptr<Transform> summary_statistics_in,
                                               std::shared_ptr<Transform> transform_function_in,
                                               const std::vector<std::string> &measurement_variables_in,
                                               const std::vector<GetMeasurementMatrixPtr> &measurement_noise_functions_in);

  virtual ~DirectGaussianMeasurementCovarianceEstimator();

  DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another);

  void operator=(const DirectGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  arma::mat get_Cygivenx() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  arma::mat get_measurement_covariance() const;
  
  void setup();
  void setup(const Parameters &parameters);

protected:
  
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);

  void make_copy(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  std::shared_ptr<Transform> transform_function;
  std::vector<GetMeasurementMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;

};

#endif
