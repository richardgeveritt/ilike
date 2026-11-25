#ifndef MEASUREMENTCOVARIANCEESTIMATOR_H
#define MEASUREMENTCOVARIANCEESTIMATOR_H

#include <memory>
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"
#include "distributions.h"
#include "ilike_header.h"

namespace ilike
{

class MeasurementCovarianceEstimator;
class MeasurementCovarianceEstimatorOutput;

// For the ensemble adjustment update (get_adjustment): V^T (I + Yhat S^-1
// Yhat^T)^-1 V, with V^T (r x J) from the compact SVD of the ensemble
// deviations and Yhat (J x d_y) the scaled deviations of the measurement
// states. By the Woodbury identity (I + Y S^-1 Y^T)^-1 = I - Y (S + Y^T Y)^-1
// Y^T, so no J x J matrix is formed (J, the ensemble size, may be large).
inline arma::mat adjustment_inner_matrix(const arma::mat &Vtranspose,
                                         const arma::mat &Yhat,
                                         const arma::mat &S)
{
  arma::mat VY = Vtranspose*Yhat;
  arma::mat output = Vtranspose*Vtranspose.t() -
    VY*arma::solve(arma::symmatu(S + Yhat.t()*Yhat), VY.t());
  return (output + output.t())/2.0;
}
class GenericMeasurementCovarianceEstimatorOutput;
class DirectGaussianMeasurementCovarianceEstimatorOutput;
class GaussianMeasurementCovarianceEstimatorOutput;
class MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
class Transform;

class MeasurementCovarianceEstimator
{
  
public:
  
  MeasurementCovarianceEstimator();
  
  // (C270) draw from this generator and seed counter instead (a copy run in
  // a parallel chunk group)
  void redirect_random_numbers(RandomNumberGenerator* rng_in, size_t* seed_in) { this->rng = rng_in; this->seed = seed_in; }
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

  // (C279) compact ensembles: whether a member's predicted measurement is a
  // function of its parameters and its noise a draw independent of them (so
  // both can be kept as rows of matrices instead of an output object per
  // member); the predicted measurement for given parameters; and the noise
  // drawn with the move (from this estimator's generator, in the order the
  // outputs' simulate draws it)
  virtual bool supports_compact_members() const { return false; }
  virtual arma::colvec compact_measurement_state(const Parameters &parameters) { Rcpp::stop("MeasurementCovarianceEstimator::compact_measurement_state - not available."); }
  virtual arma::colvec compact_noise() { Rcpp::stop("MeasurementCovarianceEstimator::compact_noise - not available."); }
  // (C281) a member's predicted measurement and noise for these parameters,
  // drawn as its output's simulate draws them (by default: the two
  // functions above, in that order); and its shift terms from them, as its
  // output's get_shift computes them (by default: the measurement plus the
  // noise times sqrt(inverse incremental temperature))
  virtual void compact_simulate(const Parameters &parameters,
                                arma::colvec &measurement_state,
                                arma::colvec &noise)
  {
    measurement_state = this->compact_measurement_state(parameters);
    noise = this->compact_noise();
  }
  virtual arma::colvec compact_shift(const arma::colvec &measurement_state,
                                     const arma::colvec &noise,
                                     double inverse_incremental_temperature) const
  {
    return measurement_state + sqrt(inverse_incremental_temperature)*noise;
  }
  // (C280) the weight a member with this predicted measurement gets at an
  // incremental temperature (a member output's evaluate_ensemble_likelihood_ratio)
  virtual double ensemble_likelihood_ratio(const arma::colvec &measurement_state,
                                           double inverse_incremental_temperature,
                                           const arma::mat &inv_sigma_precomp,
                                           double log_det_precomp) const { Rcpp::stop("MeasurementCovarianceEstimator::ensemble_likelihood_ratio - not available."); }
  
  MeasurementCovarianceEstimatorOutput* initialise();
  MeasurementCovarianceEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  virtual void setup()=0;
  virtual void setup(const Parameters &parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                              arma::mat &inv_sigma_precomp,
                                              double &log_det_precomp)=0;

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
  
  std::shared_ptr<Transform> transform;
  
  std::shared_ptr<Transform> summary_statistics;
  
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator()=0;
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)=0;
  
  virtual void setup_measurement_variables()=0;
  virtual void setup_measurement_variables(const Parameters &conditioned_on_parameters)=0;
  
  void make_copy(const MeasurementCovarianceEstimator &another);
  
  bool set_using_parameters;
  
  std::vector<std::string> measurement_variables;
  
};
}

#endif
