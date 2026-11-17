#ifndef DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTLINEARGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "direct_gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{

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
                                                     const std::vector<std::string> &state_variables_in);
  
  DirectLinearGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                     size_t* seed_in,
                                                     Data* data_in,
                                                     std::shared_ptr<Transform> transform_in,
                                                     std::shared_ptr<Transform> summary_statistics_in,
                                                     GetMatrixPtr measurement_matrix_in,
                                                     GetMatrixPtr measurement_covariance_in,
                                                     const std::string &measurement_variable_in,
                                                     const std::vector<std::string> &state_variables_in);
  
  virtual ~DirectLinearGaussianMeasurementCovarianceEstimator();
  
  DirectLinearGaussianMeasurementCovarianceEstimator(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
  
  void operator=(const DirectLinearGaussianMeasurementCovarianceEstimator &another);
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

protected:

  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);
  
  void make_copy(const DirectLinearGaussianMeasurementCovarianceEstimator &another);

  // mean of zero
  //GaussianIndependentProposalKernel kernel;
  
  std::vector<arma::mat> As;
  std::vector<GetMatrixPtr> A_functions;
  std::vector<GetMatrixPtr> measurement_noise_functions;
  
  std::vector<std::string> state_variables;
  
};
}

#endif
