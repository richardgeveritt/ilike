#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"

namespace ilike
{
class GaussianMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{
  
public:
  
  GaussianMeasurementCovarianceEstimator();
  
  GaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                         size_t* seed_in,
                                         Data* data_in,
                                         std::shared_ptr<Transform> inverse_transform_in,
                                         std::shared_ptr<Transform> summary_statistics_in,
                                         const std::vector<std::string> &measurement_variables_in);
  
  virtual ~GaussianMeasurementCovarianceEstimator();
  
  GaussianMeasurementCovarianceEstimator(const GaussianMeasurementCovarianceEstimator &another);
  
  void operator=(const GaussianMeasurementCovarianceEstimator &another);
  virtual GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const=0;
  
  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature) const;
  
  virtual arma::mat get_measurement_covariance() const=0;
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                      arma::mat &inv_sigma_precomp,
                                      double &log_det_precomp);
  
protected:
  
  void make_copy(const GaussianMeasurementCovarianceEstimator &another);
  
};
}

#endif
