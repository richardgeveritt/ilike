#include "gaussian_measurement_covariance_estimator.h"

namespace ilike
{
GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator()
:MeasurementCovarianceEstimator()
{
}

GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                               size_t* seed_in,
                                                                               Data* data_in,
                                                                               std::shared_ptr<Transform> inverse_transform_in,
                                                                               std::shared_ptr<Transform> summary_statistics_in,
                                                                               const std::vector<std::string> &measurement_variables_in)
:MeasurementCovarianceEstimator(rng_in,seed_in,data_in,inverse_transform_in,summary_statistics_in,measurement_variables_in)
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

// (C280) the ratio a member output computes (moved here so that a compact
// ensemble's rows use the same arithmetic)
double GaussianMeasurementCovarianceEstimator::ensemble_likelihood_ratio(const arma::colvec &measurement_state,
                                                                        double inverse_incremental_temperature,
                                                                        const arma::mat &inv_sigma_precomp,
                                                                        double log_det_precomp) const
{
  double c = (double(inv_sigma_precomp.n_rows)/2.0)*(1.0/inverse_incremental_temperature)*log(inverse_incremental_temperature) + (double(inv_sigma_precomp.n_rows)/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log_det_precomp;
  return c + dmvnorm_using_precomp(this->measurement,
                                   measurement_state,
                                   inv_sigma_precomp,
                                   log_det_precomp);
}

arma::mat GaussianMeasurementCovarianceEstimator::get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                                                           double inverse_incremental_temperature) const
{
  return Cyy + inverse_incremental_temperature*this->get_measurement_covariance();
}

void GaussianMeasurementCovarianceEstimator::change_data()
{
  this->current_data = this->data;
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    for (size_t i=1; i<this->measurement_variables.size(); ++i)
    {
      this->measurement = join_cols(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
    }
  }
}

void GaussianMeasurementCovarianceEstimator::change_data(std::shared_ptr<Data> new_data)
{
  this->current_data = new_data.get();
  
  if (this->measurement_variables.size()>0)
  {
    this->measurement = (*this->current_data)[this->measurement_variables[0]].as_col();
    
    if (this->measurement_variables.size()>0)
    {
      for (size_t i=1; i<this->measurement_variables.size(); ++i)
      {
        this->measurement = join_cols(this->measurement,(*this->current_data)[this->measurement_variables[i]].as_col());
      }
    }
  }
}

void GaussianMeasurementCovarianceEstimator::precompute_gaussian_covariance(double inverse_incremental_temperature,
                                                                            arma::mat &inv_sigma_precomp,
                                                                            double &log_det_precomp)
{
  arma::mat for_precomp;
  if (isinf(inverse_incremental_temperature))
  {
    arma::mat meas_cov = this->get_measurement_covariance();
    inv_sigma_precomp.zeros(meas_cov.n_rows,meas_cov.n_cols);
    log_det_precomp = arma::datum::inf;
  }
  else
  {
    for_precomp = inverse_incremental_temperature*this->get_measurement_covariance();
    inv_sigma_precomp = arma::inv_sympd(for_precomp);
    log_det_precomp = arma::log_det_sympd(for_precomp);
  }

}
}
