#include "gaussian_measurement_covariance_estimator.h"

GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator()
  :MeasurementCovarianceEstimator()
{
}

GaussianMeasurementCovarianceEstimator::~GaussianMeasurementCovarianceEstimator()
{
}

GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator(const GaussianMeasurementCovarianceEstimator &another)
  :MeasurementCovarianceEstimator(another)
{
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimator::operator=(const GaussianMeasurementCovarianceEstimator &another)
{
  if(this == &another)
    return;

  MeasurementCovarianceEstimator::operator=(another);
  this->make_copy(another);
}

void GaussianMeasurementCovarianceEstimator::make_copy(const GaussianMeasurementCovarianceEstimator &another)
{
}

bool GaussianMeasurementCovarianceEstimator::need_Cxx() const
{
  return false;
}

void GaussianMeasurementCovarianceEstimator::find_Cygivenx(const arma::mat &inv_Cxx,
                                                           const arma::mat &Cxy,
                                                           const arma::mat &Cyy)
{
  // blank on purpose
}

/*
arma::mat GaussianMeasurementCovarianceEstimator::get_kalman_gain(const arma::mat &Cxy,
                                                                  const arma::mat &Cyy,
                                                                  double inverse_incremental_temperature)
{
  return Cxy*(Cyy + inverse_temperature*this->get_measurement_covariance()).i();
}
*/
