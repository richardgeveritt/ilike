#ifndef GENERICMEASUREMENTCOVARIANCEESTIMATOR_H
#define GENERICMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

namespace ilike
{
class MeasurementCovarianceEstimatorOutput;
class GenericMeasurementCovarianceEstimatorOutput;

class GenericMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{
  
public:
  
  GenericMeasurementCovarianceEstimator();
  
  GenericMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                        size_t* seed_in,
                                        Data* data_in,
                                        std::shared_ptr<Transform> transform_in,
                                        std::shared_ptr<Transform> summary_statistics_in,
                                        SimulateModelPtr simulator_in,
                                        const std::vector<std::string> &measurement_variables_in);
  
  virtual ~GenericMeasurementCovarianceEstimator();
  
  GenericMeasurementCovarianceEstimator(const GenericMeasurementCovarianceEstimator &another);
  
  void operator=(const GenericMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  
  void setup();
  void setup(const Parameters &parameters);
  
  //void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  bool need_Cxx() const;
  
  void find_Cygivenx(const arma::mat &inv_Cxx,
                     const arma::mat &Cxy,
                     const arma::mat &Cyy);
  
  arma::mat get_adjustment(const arma::mat &Zf,
                           const arma::mat &Dhathalf,
                           const arma::mat &P,
                           const arma::mat &Vtranspose,
                           const arma::mat &Yhat,
                           double inverse_incremental_temperature) const;
  
  arma::mat get_sqrt_adjustment(const arma::mat &Sigma,
                                const arma::mat &HSigmaHt,
                                double inverse_incremental_temperature) const;
  
  arma::mat get_Cygivenx() const;
  
  arma::mat get_unconditional_measurement_covariance(const arma::mat &Cyy,
                                                     double inverse_incremental_temperature) const;
  
  Parameters simulate(const Parameters &current_state);
  
  arma::colvec gaussian_simulate();
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                      arma::mat &inv_sigma_precomp,
                                      double &log_det_precomp);
  
protected:
  
  friend GenericMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);
  
  void make_copy(const GenericMeasurementCovarianceEstimator &another);
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  //GetMeasurementMatrixPtr measurement_noise_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //arma::mat measurement_noise;
  
  SimulateModelPtr simulator;
  
  GaussianIndependentProposalKernel gaussian_simulator;
  
  size_t dimension;
  
  arma::mat Cygivenx;
  
};
}

#endif
