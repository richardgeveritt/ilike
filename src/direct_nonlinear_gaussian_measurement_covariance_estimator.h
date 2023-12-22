#ifndef DIRECTNONLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTNONLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

//class DirectGaussianMeasurementCovarianceEstimatorOutput;

class DirectNonLinearGaussianMeasurementCovarianceEstimator : public DirectGaussianMeasurementCovarianceEstimator
{

public:

  DirectNonLinearGaussianMeasurementCovarianceEstimator();

  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const arma::mat &measurement_noise_in);
                                                        //const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<arma::mat> &measurement_noise_in);
                                                        //const std::vector<std::string> &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        GetMatrixPtr measurement_noise_function_in);
                                                        //const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<GetMatrixPtr> &measurement_noise_functions_in);
                                                        //const std::vector<std::string> &measurements_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const arma::mat &measurement_noise_in,
  const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<arma::mat> &measurement_noise_in,
  const std::vector<std::string> &measurement_variables_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        GetMatrixPtr measurement_noise_function_in,
  const std::string &measurement_variable_in);
  
  DirectNonLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                        size_t* seed_in,
                                                        Data* data_in,
                                                        std::shared_ptr<Transform> transform_in,
                                                        std::shared_ptr<Transform> summary_statistics_in,
                                                        std::shared_ptr<Transform> transform_function_in,
                                                        const std::vector<GetMatrixPtr> &measurement_noise_functions_in,
  const std::vector<std::string> &measurements_variables_in);

  virtual ~DirectNonLinearGaussianMeasurementCovarianceEstimator();

  DirectNonLinearGaussianMeasurementCovarianceEstimator(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);

  void operator=(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);
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

  void make_copy(const DirectNonLinearGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  std::shared_ptr<Transform> transform_function;
  std::vector<GetMatrixPtr> measurement_noise_functions;
  //arma::mat measurement_noise;

};

#endif
