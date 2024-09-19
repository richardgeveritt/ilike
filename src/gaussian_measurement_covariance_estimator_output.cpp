#include "gaussian_measurement_covariance_estimator_output.h"
#include "gaussian_measurement_covariance_estimator.h"
#include "distributions.h"

namespace ilike
{
GaussianMeasurementCovarianceEstimatorOutput::GaussianMeasurementCovarianceEstimatorOutput()
:MeasurementCovarianceEstimatorOutput()
{
}

GaussianMeasurementCovarianceEstimatorOutput::~GaussianMeasurementCovarianceEstimatorOutput()
{
}

GaussianMeasurementCovarianceEstimatorOutput::GaussianMeasurementCovarianceEstimatorOutput(const GaussianMeasurementCovarianceEstimatorOutput &another)
:MeasurementCovarianceEstimatorOutput(another)
{
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimatorOutput::operator=(const GaussianMeasurementCovarianceEstimatorOutput &another)
{
  if(this == &another)
    return;
  
  MeasurementCovarianceEstimatorOutput::operator=(another);
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimatorOutput::make_copy(const GaussianMeasurementCovarianceEstimatorOutput &another)
{
  this->measurement_state = another.measurement_state;
  this->random_shift = another.random_shift;
}

arma::rowvec GaussianMeasurementCovarianceEstimatorOutput::get_measurement_state_for_covariance() const
{
  return arma::conv_to<arma::rowvec>::from(this->measurement_state);
}

arma::colvec GaussianMeasurementCovarianceEstimatorOutput::get_shift(double inverse_incremental_temperature) const
{
  return this->measurement_state + sqrt(inverse_incremental_temperature)*this->random_shift;
}

arma::colvec GaussianMeasurementCovarianceEstimatorOutput::get_deterministic_shift() const
{
  return this->measurement_state;
}

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                        const arma::mat &inv_sigma_precomp,
                                                                                        double log_det_precomp) const
{
  return dmvnorm_using_precomp(*this->get_gaussian_estimator()->get_measurement_pointer(),
                               this->measurement_state,
                               inv_sigma_precomp,
                               log_det_precomp);
}

/*
 double GaussianMeasurementCovarianceEstimatorOutput::evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 inverse_incremental_temperature*this->get_gaussian_estimator()->get_measurement_covariance());
 }
 */

double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                                                                  const arma::mat &inv_sigma_precomp,
                                                                                                  double log_det_precomp) const
{
  // parameters of covariance should already be set at this point, so second argument does nothing
  //return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
  //               this->measurement_state,
  //               inverse_incremental_temperature*this->get_gaussian_estimator()->get_measurement_covariance());
  
  return dmvnorm_using_precomp(*this->get_gaussian_estimator()->get_measurement_pointer(),
                               this->measurement_state,
                               inv_sigma_precomp,
                               log_det_precomp);
}

/*
 double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
 const Parameters &conditioned_on_parameters)
 {
 // parameters of covariance should already be set at this point, so second argument does nothing
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 inverse_incremental_temperature*this->get_gaussian_estimator()->get_measurement_covariance());
 }
 */

double GaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_gaussian_estimator()->get_measurement_covariance());
}

/*
 double GaussianMeasurementCovarianceEstimatorOutput::evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 this->get_gaussian_estimator()->get_measurement_covariance());
 }
 */

double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood()
{
  return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
                 this->measurement_state,
                 this->get_gaussian_estimator()->get_measurement_covariance());
}

/*
 double GaussianMeasurementCovarianceEstimatorOutput::subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)
 {
 return dmvnorm(*this->get_estimator()->get_measurement_pointer(),
 this->measurement_state,
 this->get_gaussian_estimator()->get_measurement_covariance());
 }
 */

void GaussianMeasurementCovarianceEstimatorOutput::write_to_file(const std::string &directory_name,
                                                                 const std::string &index)
{
  
}

void GaussianMeasurementCovarianceEstimatorOutput::close_ofstreams()
{
  
}

}
