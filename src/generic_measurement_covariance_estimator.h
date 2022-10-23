#ifndef GENERICMEASUREMENTCOVARIANCEESTIMATOR_H
#define GENERICMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"
#include "function_pointers.h"
#include "gaussian_independent_proposal_kernel.h"

class MeasurementCovarianceEstimatorOutput;
class GenericMeasurementCovarianceEstimatorOutput;

class GenericMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{

public:

  GenericMeasurementCovarianceEstimator();
  virtual ~GenericMeasurementCovarianceEstimator();

  GenericMeasurementCovarianceEstimator(const GenericMeasurementCovarianceEstimator &another);

  void operator=(const GenericMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
  arma::mat get_Cygivenx() const;
  
  arma::colvec simulate(const Parameters &current_state);
  
  arma::colvec gaussian_simulate();

protected:
  
  friend GenericMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);

  void make_copy(const GenericMeasurementCovarianceEstimator &another);
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  //GetMeasurementMatrixPtr measurement_noise_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //arma::mat measurement_noise;
  
  SimulateIndependentProposalPtr simulator;
  
  GaussianIndependentProposalKernel gaussian_simulator;
  
  size_t dimension;
  
  arma::mat Cygivenx;
  
};

#endif
