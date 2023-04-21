#ifndef MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H
#define MEASUREMENTCOVARIANCEESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"

class Parameters;
class MeasurementCovarianceEstimator;
class Index;

class MeasurementCovarianceEstimatorOutput
{

public:

  MeasurementCovarianceEstimatorOutput();
  virtual ~MeasurementCovarianceEstimatorOutput();

  MeasurementCovarianceEstimatorOutput(const MeasurementCovarianceEstimatorOutput &another);

  void operator=(const MeasurementCovarianceEstimatorOutput &another);
  virtual MeasurementCovarianceEstimatorOutput* duplicate() const=0;
  
  void simulate(const Parameters &parameters);
  void subsample_simulate(const Parameters &parameters);
  virtual void specific_simulate(const Parameters &parameters)=0;
  virtual void subsample_specific_simulate(const Parameters &parameters)=0;
  virtual arma::rowvec get_measurement_state_for_covariance() const=0;
  virtual arma::colvec get_shift(double inverse_incremental_temperature) const=0;
  virtual arma::colvec get_deterministic_shift() const=0;
  
  arma::colvec* get_measurement();
  
  virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature)=0;
  /*
  virtual double evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                    const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature)=0;
  /*
  virtual double subsample_evaluate_ensemble_likelihood_ratio(double inverse_incremental_temperature,
                                                              const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual double evaluate_likelihood()=0;
  //virtual double evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_likelihood()=0;
  //virtual double subsample_evaluate_likelihood(const Parameters &conditioned_on_parameters)=0;
  
  virtual MeasurementCovarianceEstimator* get_estimator()=0;
  
  virtual void close_ofstreams()=0;
  
  void write(const std::string &directory_name);
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index="")=0;
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
protected:
  
  bool write_to_file_flag;
  
  // Stored in ModelAndAlgorithm or in main.
  //MeasurementCovarianceEstimator* estimator;
  
  void make_copy(const MeasurementCovarianceEstimatorOutput &another);

};

#endif
