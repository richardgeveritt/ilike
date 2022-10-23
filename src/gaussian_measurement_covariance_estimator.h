#ifndef GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define GAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"

class GaussianMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{

public:

  GaussianMeasurementCovarianceEstimator();
  virtual ~GaussianMeasurementCovarianceEstimator();

  GaussianMeasurementCovarianceEstimator(const GaussianMeasurementCovarianceEstimator &another);

  void operator=(const GaussianMeasurementCovarianceEstimator &another);
  virtual GaussianMeasurementCovarianceEstimator* gaussian_duplicate() const=0;
  
  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
  //virtual arma::mat get_measurement_covariance() const=0;

protected:

  void make_copy(const GaussianMeasurementCovarianceEstimator &another);
  
};

#endif
