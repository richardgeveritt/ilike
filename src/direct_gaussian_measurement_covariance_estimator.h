#ifndef DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define DIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "gaussian_measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{

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
                                               const std::vector<std::string> &measurement_variables_in);
  
  virtual ~DirectGaussianMeasurementCovarianceEstimator();
  
  DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  void operator=(const DirectGaussianMeasurementCovarianceEstimator &another);
  
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  // (C264: virtual, so the linear case can stack its values without a map)
  virtual arma::colvec get_measurement_state(const Parameters &parameters) const;
  virtual Parameters get_measurement_state_parameters(const Parameters &parameters) const=0;

  arma::mat get_Cygivenx() const;

  // (C279) see MeasurementCovarianceEstimator; the output's simulate uses the
  // same two functions
  bool supports_compact_members() const { return this->transform == NULL; }
  arma::colvec compact_measurement_state(const Parameters &parameters);
  arma::colvec compact_noise();

  arma::mat get_measurement_covariance() const;
  
  void setup();
  void setup(const Parameters &parameters);
  
protected:
  
  // mean of zero
  GaussianIndependentProposalKernel kernel;
  
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);

  void make_copy(const DirectGaussianMeasurementCovarianceEstimator &another);

  // mean of zero
  //GaussianIndependentProposalKernel kernel;

};
}

#endif
