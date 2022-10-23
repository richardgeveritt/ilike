#ifndef MEASUREMENTCOVARIANCEESTIMATOR_H
#define MEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"
#include "distributions.h"

class MeasurementCovarianceEstimator;
class MeasurementCovarianceEstimatorOutput;
class Data;

class MeasurementCovarianceEstimator
{

public:

  MeasurementCovarianceEstimator();
  virtual ~MeasurementCovarianceEstimator();

  MeasurementCovarianceEstimator(const MeasurementCovarianceEstimator &another);

  void operator=(const MeasurementCovarianceEstimator &another);
  virtual MeasurementCovarianceEstimator* duplicate() const=0;
  
  MeasurementCovarianceEstimatorOutput* initialise();
  MeasurementCovarianceEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  virtual bool need_Cxx() const=0;
  
  virtual void find_Cygivenx(const arma::mat &inv_Cxx,
                             const arma::mat &Cxy,
                             const arma::mat &Cyy)=0;
  
  void change_data();
  void change_data(Data* new_data);
  
  Data* get_data() const;
  
  arma::colvec* get_measurement_pointer();

protected:
  
  // Not stored here. Stored in "main'.
  Data* data;
  
  // not stored here
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  size_t* seed;
  
  arma::colvec measurement;
  
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator()=0;
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)=0;
  
  void make_copy(const MeasurementCovarianceEstimator &another);
  
  bool set_using_parameters;
  
  std::vector<std::string> measurement_variables;

};

#endif
