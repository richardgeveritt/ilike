#include "direct_gaussian_measurement_covariance_estimator.h"
#include "direct_gaussian_measurement_covariance_estimator_output.h"
#include "transform.h"

namespace ilike
{
DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator()
:GaussianMeasurementCovarianceEstimator()
{
}

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                                           size_t* seed_in,
                                                                                           Data* data_in,
                                                                                           std::shared_ptr<Transform> transform_in,
                                                                                           std::shared_ptr<Transform> summary_statistics_in,
                                                                                           const std::vector<std::string> &measurement_variables_in)
: GaussianMeasurementCovarianceEstimator(rng_in,
                                         seed_in,
                                         data_in,
                                         transform_in,
                                         summary_statistics_in,
                                         measurement_variables_in)
{
  
}

DirectGaussianMeasurementCovarianceEstimator::~DirectGaussianMeasurementCovarianceEstimator()
{
}

DirectGaussianMeasurementCovarianceEstimator::DirectGaussianMeasurementCovarianceEstimator(const DirectGaussianMeasurementCovarianceEstimator &another)
:GaussianMeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void DirectGaussianMeasurementCovarianceEstimator::operator=(const DirectGaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;
  
  GaussianMeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

void DirectGaussianMeasurementCovarianceEstimator::make_copy(const DirectGaussianMeasurementCovarianceEstimator &another)
{
  this->kernel = another.kernel;
}

MeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new DirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* DirectGaussianMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new DirectGaussianMeasurementCovarianceEstimatorOutput(this);
  return output;
}

void DirectGaussianMeasurementCovarianceEstimator::setup()
{
  this->setup_measurement_variables();
}

void DirectGaussianMeasurementCovarianceEstimator::setup(const Parameters &parameters)
{
  this->setup_measurement_variables(parameters);
}

arma::mat DirectGaussianMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::mat DirectGaussianMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->kernel.get_covariance(this->measurement_variables);
}

arma::colvec DirectGaussianMeasurementCovarianceEstimator::get_measurement_state(const Parameters &parameters) const
{
  return this->get_measurement_state_parameters(parameters).get_colvec(this->measurement_variables);
}
}
