#include "gaussian_measurement_covariance_estimator.h"

GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator()
  :MeasurementCovarianceEstimator()
{
}

GaussianMeasurementCovarianceEstimator::GaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                                               size_t* seed_in,
                                                                               Data* data_in,
                                                                               std::shared_ptr<Transform> inverse_transform_in,
                                                                               std::shared_ptr<Transform> summary_statistics_in)
:MeasurementCovarianceEstimator(rng_in,seed_in,data_in,inverse_transform_in,summary_statistics_in)
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

arma::mat GaussianMeasurementCovarianceEstimator::get_adjustment(const arma::mat &Zf,
                                                                 const arma::mat &Ginv,
                                                                 const arma::mat &Ftranspose,
                                                                 const arma::mat &V,
                                                                 double inverse_incremental_temperature)
{
  // follows https://arxiv.org/abs/2006.02941
  arma::mat for_eig = V*(inverse_incremental_temperature*this->get_measurement_covariance())*V.t();
  
  arma::mat C;
  arma::vec diagGamma;
  arma::mat Ctrans;
  arma::svd(C,diagGamma,Ctrans,for_eig);
  
  arma::mat Gamma(diagGamma.n_elem,diagGamma.n_elem);
  Gamma.diag() = diagGamma;
  arma::mat I;
  I.eye(diagGamma.n_elem,diagGamma.n_elem);
  
  return Zf*C*arma::sqrtmat_sympd(arma::inv_sympd(I+Gamma))*Ginv*Ftranspose;
}

arma::mat GaussianMeasurementCovarianceEstimator::get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                                                           double inverse_incremental_temperature)
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

void GaussianMeasurementCovarianceEstimator::change_data(Data* new_data)
{
  this->current_data = new_data;
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

void GaussianMeasurementCovarianceEstimator::precompute_gaussian_covariance(double inverse_incremental_temperature)
{
  arma::mat for_precomp = inverse_incremental_temperature*this->get_measurement_covariance();
  this->inv_sigma_precomp = arma::inv_sympd(for_precomp);
  this->log_det_precomp = arma::log_det_sympd(for_precomp);
}
