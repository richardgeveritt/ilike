#ifndef DIRECTABCGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTABCGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

//class DirectGaussianMeasurementCovarianceEstimatorOutput;

namespace ilike
{
class DirectABCGaussianMeasurementCovarianceEstimator : public DirectGaussianMeasurementCovarianceEstimator
{
  
public:
  
  DirectABCGaussianMeasurementCovarianceEstimator();
  
  
  DirectABCGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                  size_t* seed_in,
                                                  Data* data_in,
                                                  std::shared_ptr<Transform> transform_in,
                                                  std::shared_ptr<Transform> summary_statistics_in,
                                                  double min_epsilon_in,
                                                  //const std::string &tempering_variable_in,
                                                  const std::string &scale_variable_in,
                                                  const std::vector<std::string> &measurement_variables_in);
  
  virtual ~DirectABCGaussianMeasurementCovarianceEstimator();
  
  DirectABCGaussianMeasurementCovarianceEstimator(const DirectABCGaussianMeasurementCovarianceEstimator &another);
  
  void operator=(const DirectABCGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const;
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  Parameters get_measurement_state_parameters(const Parameters &parameters) const;
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Dhathalf,
                           const arma::mat &P,
                           const arma::mat &Vtranspose,
                           const arma::mat &Yhat,
                           double inverse_incremental_temperature) const;
  
  arma::mat get_sqrt_adjustment(const arma::mat &Sigma,
                                const arma::mat &HSigmaHt,
                                double inverse_incremental_temperature) const;
  
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
  
  void make_copy(const DirectABCGaussianMeasurementCovarianceEstimator &another);
  
  //Parameters conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //std::shared_ptr<Transform> summary_statistics;
  double min_epsilon;
  //arma::mat measurement_noise;
  
  //std::string tempering_variable;
  std::string scale_variable;
  
};
}

#endif
