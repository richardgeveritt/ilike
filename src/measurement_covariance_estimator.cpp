#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "data.h"

MeasurementCovarianceEstimator::MeasurementCovarianceEstimator()
{
  // when we read in the data, get the vector corresponding to the measurement name and vectorise it, and store it in
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
  this->measurement = another.measurement;
}

MeasurementCovarianceEstimatorOutput* MeasurementCovarianceEstimator::initialise()
{
  return this->initialise_measurement_covariance_estimator();
}

MeasurementCovarianceEstimatorOutput* MeasurementCovarianceEstimator::initialise(const Parameters &conditioned_on_parameters)
{
  return this->initialise_measurement_covariance_estimator(conditioned_on_parameters);
}

void MeasurementCovarianceEstimator::change_data()
{
  this->current_data = this->data;
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    if (this->measurement_variables.size()>0)
    {
      for (size_t i=1; i<this->measurement_variables.size(); ++i)
      {
        this->measurement = join_rows(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
      }
    }
  }
}

void MeasurementCovarianceEstimator::change_data(Data* new_data)
{
  this->current_data = new_data;
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    if (this->measurement_variables.size()>0)
    {
      for (size_t i=1; i<this->measurement_variables.size(); ++i)
      {
        this->measurement = join_rows(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
      }
    }
  }
}

Data* MeasurementCovarianceEstimator::get_data() const
{
  return this->data;
}

arma::colvec* MeasurementCovarianceEstimator::get_measurement_pointer()
{
  return &this->measurement;
}
