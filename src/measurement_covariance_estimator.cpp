#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "parameters.h"

MeasurementCovarianceEstimator::MeasurementCovarianceEstimator()
{
  this->transform = NULL;
  this->summary_statistics = NULL;
}

MeasurementCovarianceEstimator::MeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                               size_t* seed_in,
                                                               Data* data_in,
                                                               std::shared_ptr<Transform> transform_in,
                                                               std::shared_ptr<Transform> summary_statistics_in,
                                                               const std::vector<std::string> &measurement_variables_in)
{
  this->data = data_in;
  this->current_data = this->data;
  this->rng = rng_in;
  this->seed = seed_in;
  this->measurement_variables = measurement_variables_in;//this->data->get_vector_variables(); // automatically get measurement variables - might want to remove this
  this->set_using_parameters = false;
  this->transform = transform_in;
  this->summary_statistics = summary_statistics_in;
}

MeasurementCovarianceEstimator::~MeasurementCovarianceEstimator()
{
}

MeasurementCovarianceEstimator::MeasurementCovarianceEstimator(const MeasurementCovarianceEstimator &another)
{
  this->make_copy(another);
}

void MeasurementCovarianceEstimator::operator=(const MeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void MeasurementCovarianceEstimator::make_copy(const MeasurementCovarianceEstimator &another)
{
  this->set_using_parameters = another.set_using_parameters;
  this->measurement_variables = another.measurement_variables;
  this->data = another.data;
  this->current_data = another.current_data;
  this->seed = another.seed;
  this->rng = another.rng;
  this->measurement = another.measurement;
  this->transform = another.transform;
  this->summary_statistics = another.summary_statistics;
  //this->inv_sigma_precomp = another.inv_sigma_precomp;
  //this->log_det_precomp = another.log_det_precomp;
}

MeasurementCovarianceEstimatorOutput* MeasurementCovarianceEstimator::initialise()
{
  //this->setup_measurement_variables();
  return this->initialise_measurement_covariance_estimator();
}

MeasurementCovarianceEstimatorOutput* MeasurementCovarianceEstimator::initialise(const Parameters &conditioned_on_parameters)
{
  //this->setup_measurement_variables(conditioned_on_parameters);
  return this->initialise_measurement_covariance_estimator(conditioned_on_parameters);
}

Data* MeasurementCovarianceEstimator::get_data() const
{
  return this->data;
}

arma::colvec* MeasurementCovarianceEstimator::get_measurement_pointer()
{
  return &this->measurement;
}

const arma::colvec* MeasurementCovarianceEstimator::get_measurement_pointer() const
{
  return &this->measurement;
}

std::vector<std::string> MeasurementCovarianceEstimator::get_measurement_variables() const
{
  return this->measurement_variables;
}
