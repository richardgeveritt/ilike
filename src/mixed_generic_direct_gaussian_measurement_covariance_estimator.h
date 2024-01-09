#ifndef MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H
#define MIXEDGENERICDIRECTGAUSSIANMEASUREMENTCOVARIANCEESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "measurement_covariance_estimator.h"
#include "ilike_header.h"
#include "gaussian_independent_proposal_kernel.h"

class MeasurementCovarianceEstimatorOutput;
class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;

class MixedGenericDirectGaussianMeasurementCovarianceEstimator : public MeasurementCovarianceEstimator
{

public:

  MixedGenericDirectGaussianMeasurementCovarianceEstimator();
  
  MixedGenericDirectGaussianMeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                                           size_t* seed_in,
                                                           Data* data_in,
                                                           std::shared_ptr<Transform> transform_in,
                                                           std::shared_ptr<Transform> summary_statistics_in,
                                                           Data* prior_data_in,
                                                           SimulateModelPtr simulator_in,
                                                           const std::vector<std::string> &prior_measurement_variables_in,
                                                           const std::vector<arma::mat> &prior_measurement_noises_in);
  
  virtual ~MixedGenericDirectGaussianMeasurementCovarianceEstimator();

  MixedGenericDirectGaussianMeasurementCovarianceEstimator(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);

  void operator=(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);
  MeasurementCovarianceEstimator* duplicate() const;
  
  void setup();
  void setup(const Parameters &parameters);
  
  void set_parameters(const Parameters &conditioned_on_parameters_in);
  
  //arma::mat get_measurement_covariance() const;
  
  bool need_Cxx() const;
  
  void find_partial_Cygivenx(const arma::mat &Cxy,
                             const arma::mat &Cyy,
                             const arma::mat &packed_members);
  
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
  
  arma::mat get_prior_measurement_covariance() const;
  arma::mat get_prior_measurement_covariance_embedded_in_full_space() const;
  arma::mat get_measurement_covariance_for_likelihood_ratio(double inverse_incremental_temperature) const;
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  
  void precompute_gaussian_covariance(double inverse_incremental_temperature);

protected:
  
  friend MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
  
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator();
  MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters);
  
  void setup_measurement_variables();
  void setup_measurement_variables(const Parameters &conditioned_on_parameters);

  void make_copy(const MixedGenericDirectGaussianMeasurementCovarianceEstimator &another);
  
  // not stored here
  //const Parameters* conditioned_on_parameters;
  
  //GetSimulateMeasurementKernelPtr measurement_kernel_function;
  //GetMeasurementMatrixPtr measurement_noise_function;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  //arma::mat measurement_noise;
  
  SimulateModelPtr simulator;
  
  GaussianIndependentProposalKernel gaussian_simulator;
  
  size_t measurement_dimension; // dimension of prior
  
  arma::mat Cygivenx;
  arma::mat likelihood_Cygivenx;
  
  // mean of zero
  GaussianIndependentProposalKernel kernel;
  
  //SimulateMeasurementKernelPtr measurement_kernel;
  TransformToDataPtr transform_function;
  std::vector<GetMatrixPtr> prior_measurement_noise_functions;
  
  // in this class, this->measurement_variables (in base) always refers to the "true" measurement variables, for the likelihood term - not those for the prior, which we store in this derived class (same goes for the prior "Data")
  std::vector<std::string> prior_measurement_variables;
  
  // not stored here - points to prior state mean in base
  Data* prior_data;
  
};

#endif
