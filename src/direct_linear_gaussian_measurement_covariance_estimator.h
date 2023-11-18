#ifndef DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

//class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectLinearGaussianMeasurementCovarianceEstimator : public DirectGaussianMeasurementCovarianceEstimator
{

public:

  DirectLinearGaussianMeasurementCovarianceEstimator();
  
  DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                     size_t* seed_in,
                                                     Data* data_in,
                                                     std::shared_ptr<Transform> transform_in,
                                                     std::shared_ptr<Transform> summary_statistics_in,
                                                     const arma::mat &measurement_matrix_in,
                                                     const arma::mat &measurement_covariance_in,
                                                     const std::string &measurement_variable_in,
                                                     const std::string &state_variable_in);
  
  DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                     size_t* seed_in,
                                                     Data* data_in,
                                                     std::shared_ptr<Transform> transform_in,
                                                     std::shared_ptr<Transform> summary_statistics_in,
                                                     GetMatrixPtr measurement_matrix_in,
                                                     GetMatrixPtr measurement_covariance_in,
                                                     const std::string &measurement_variable_in,
                                                     const std::string &state_variable_in);

  virtual ~DirectLinearGaussianMeasurementCovarianceEstimator();

  DirectLinearGaussianMeasurementCovarianceEstimator(const DirectLinearGaussianMeasurementCovarianceEstimator &another);

  void operator=(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  Parameters get_measurement_state_parameters(const Parameters &parameters) const;
  
  //GaussianIndependentProposalKernel get_kernel() const;
  
  //arma::mat get_Cygivenx() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  //void setup();
  //void setup(const Parameters &parameters);

protected:
  
  //friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
  //MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  //MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);

  void make_copy(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  std::vector<arma::mat> As;
  std::vector<GetMatrixPtr> A_functions;
  std::vector<GetMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;
  
  std::string state_variable;

};

#endif
