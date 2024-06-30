#ifndef MEASUREMENTCOVARIANCEESTIMATOR_H
#define MEASUREMENTCOVARIANCEESTIMATOR_H

#include <memory>
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"
#include "distributions.h"
#include "ilike_header.h"

class MeasurementCovarianceEstimator;
class MeasurementCovarianceEstimatorOutput;
class GenericMeasurementCovarianceEstimatorOutput;
class DirectGaussianMeasurementCovarianceEstimatorOutput;
class GaussianMeasurementCovarianceEstimatorOutput;
class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
class Transform;

class MeasurementCovarianceEstimator
{

public:

  MeasurementCovarianceEstimator();
  
  MeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                 size_t* seed_in,
                                 Data* data_in,
                                 std::shared_ptr<Transform> transform_in,
                                 std::shared_ptr<Transform> summary_statistics_in,
                                 const std::vector<std::string> &measurement_variables_in);
  
  virtual ~MeasurementCovarianceEstimator();

  MeasurementCovarianceEstimator(const MeasurementCovarianceEstimator &another);

  void operator=(const MeasurementCovarianceEstimator &another);
  virtual MeasurementCovarianceEstimator* duplicate() const=0;
  
  MeasurementCovarianceEstimatorOutput* initialise();
  MeasurementCovarianceEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  virtual void setup()=0;
  virtual void setup(const Parameters &parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                              arma::mat &inv_sigma_precomp,
                                              double &log_det_precomp)=0;
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const arma::mat &Cxy,
                             const arma::mat &Cyy)=0;
  
  virtual arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                             double inverse_incremental_temperature) const=0;
  
  virtual arma::mat get_adjustment(const arma::mat &Zf,
                                   const arma::mat &Dhathalf,
                                   const arma::mat &P,
                                   const arma::mat &Vtranspose,
                                   const arma::mat &Yhat,
                                   double inverse_incremental_temperature) const=0;
  
  virtual arma::mat get_sqrt_adjustment(const arma::mat &Sigma,
                                        const arma::mat &HSigmaHt,
                                        double inverse_incremental_temperature) const=0;
  
  virtual arma::mat get_Cygivenx() const=0;
  
  virtual void change_data()=0;
  virtual void change_data(std::shared_ptr<Data> new_data)=0;
  
  Data* get_data() const;
  
  arma::colvec* get_measurement_pointer();
  const arma::colvec* get_measurement_pointer() const;
  
  std::vector<std::string> get_measurement_variables() const;

protected:
  
  friend MeasurementCovarianceEstimatorOutput;
  friend GenericMeasurementCovarianceEstimatorOutput;
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  friend GaussianMeasurementCovarianceEstimatorOutput;
  friend MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
  
  // Not stored here. Stored in "main'.
  Data* data;
  
  // not stored here
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  size_t* seed;
  
  arma::colvec measurement;
  
  //TransformPtr inverse_transform;
  std::shared_ptr<Transform> transform;
  
  //SummaryStatisticsPtr summary_statistics;
  std::shared_ptr<Transform> summary_statistics;
  
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator()=0;
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)=0;
  
  virtual void setup_measurement_variables()=0;
  virtual void setup_measurement_variables(const Parameters &conditioned_on_parameters)=0;
  
  void make_copy(const MeasurementCovarianceEstimator &another);
  
  bool set_using_parameters;
  
  std::vector<std::string> measurement_variables;

};

#endif
