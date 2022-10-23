#include "generic_measurement_covariance_estimator.h"
#include "generic_measurement_covariance_estimator_output.h"

GenericMeasurementCovarianceEstimator::GenericMeasurementCovarianceEstimator()
  :MeasurementCovarianceEstimator()
{
}

GenericMeasurementCovarianceEstimator::~GenericMeasurementCovarianceEstimator()
{
}

GenericMeasurementCovarianceEstimator::GenericMeasurementCovarianceEstimator(const GenericMeasurementCovarianceEstimator &another)
  :MeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimator::operator=(const GenericMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  MeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

void GenericMeasurementCovarianceEstimator::make_copy(const GenericMeasurementCovarianceEstimator &another)
{
  //this->conditioned_on_parameters = another.conditioned_on_parameters;
  //this->measurement_kernel = another.measurement_kernel;
  //this->measurement_noise = another.measurement_noise;
  //this->measurement_noise_function = another.measurement_noise_function;
  //this->measurement_kernel_function = another.measurement_kernel_function;
  this->dimension = another.dimension;
  this->simulator = another.simulator;
  this->Cygivenx = another.Cygivenx;
  this->gaussian_simulator = another.gaussian_simulator;
}

MeasurementCovarianceEstimator* GenericMeasurementCovarianceEstimator::duplicate() const
{
  return( new GenericMeasurementCovarianceEstimator(*this));
}


MeasurementCovarianceEstimatorOutput* GenericMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator()
{
  MeasurementCovarianceEstimatorOutput* output = new GenericMeasurementCovarianceEstimatorOutput(this);
  return output;
}

MeasurementCovarianceEstimatorOutput* GenericMeasurementCovarianceEstimator::initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)
{
  MeasurementCovarianceEstimatorOutput* output = new GenericMeasurementCovarianceEstimatorOutput(this);
  return output;
}

/*
arma::mat GenericMeasurementCovarianceEstimator::get_measurement_covariance() const
{
  arma::mat measurement_noise;
  measurement_noise.zeros(this->dimension,this->dimension);
  return measurement_noise;
}
*/

arma::mat GenericMeasurementCovarianceEstimator::get_Cygivenx() const
{
  return this->Cygivenx;
}

/*
void GenericMeasurementCovarianceEstimator::set_parameters(const Parameters &conditioned_on_parameters_in)
{
  if (this->set_using_parameters)
  {
    this->conditioned_on_parameters = &conditioned_on_parameters_in;
    this->measurement_noise = this->measurement_noise_function(conditioned_on_parameters_in);
  }
}
 */

bool GenericMeasurementCovarianceEstimator::need_Cxx() const
{
  return true;
}

void GenericMeasurementCovarianceEstimator::find_Cygivenx(const arma::mat &inv_Cxx,
                                                          const arma::mat &Cxy,
                                                          const arma::mat &Cyy)
{
  this->Cygivenx = Cyy - Cxy.t()*inv_Cxx*Cxy;
}

arma::colvec GenericMeasurementCovarianceEstimator::simulate(const Parameters &current_state)
{
  return this->simulator(*this->rng,
                         current_state).get_vector(this->measurement_variables);
}

arma::colvec GenericMeasurementCovarianceEstimator::gaussian_simulate()
{
  return this->gaussian_simulator.independent_simulate(*this->rng).get_vector(this->measurement_variables);
}
