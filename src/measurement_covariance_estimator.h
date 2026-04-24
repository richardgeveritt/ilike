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
  /**
   * @file measurement_covariance_estimator.h
   * @brief Defines the MeasurementCovarianceEstimator class.
   *
   * Estimates the measurement covariance for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class MeasurementCovarianceEstimator
   * @brief The measurement covariance estimator class.
   */


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
  
  /**
   * @brief Default constructor for MeasurementCovarianceEstimator.
   */
  MeasurementCovarianceEstimator();
  
  MeasurementCovarianceEstimator(RandomNumberGenerator* rng_in,
                                 size_t* seed_in,
                                 Data* data_in,
                                 std::shared_ptr<Transform> transform_in,
                                 std::shared_ptr<Transform> summary_statistics_in,
                                 const std::vector<std::string> &measurement_variables_in);
  
  /**
   * @brief Destructor for MeasurementCovarianceEstimator.
   */
  virtual ~MeasurementCovarianceEstimator();
  
  /**
   * @brief Copy constructor for MeasurementCovarianceEstimator.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  MeasurementCovarianceEstimator(const MeasurementCovarianceEstimator &another);
  
  /**
   * @brief Assignment operator for MeasurementCovarianceEstimator.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  void operator=(const MeasurementCovarianceEstimator &another);
  /**
   * @brief Creates a deep copy of this MeasurementCovarianceEstimator object.
   *
   * @return The result.
   */
  virtual MeasurementCovarianceEstimator* duplicate() const=0;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* initialise();
  /**
   * @brief Initialises the estimator with given parameters and returns an output object.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  MeasurementCovarianceEstimatorOutput* initialise(const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs any setup required before running the algorithm.
   */
  virtual void setup()=0;
  /**
   * @brief Performs setup given the supplied parameters.
   *
   * @param parameters The parameters.
   */
  virtual void setup(const Parameters &parameters)=0;
  
  virtual void precompute_gaussian_covariance(double inverse_incremental_temperature,
                                              arma::mat &inv_sigma_precomp,
                                              double &log_det_precomp)=0;
  
  //virtual arma::mat get_measurement_covariance()=0;
  //virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
  /**
   * @brief Performs the need cxx operation.
   *
   * @return The result.
   */
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
  
  /**
   * @brief Returns the cygivenx.
   *
   * @return The result.
   */
  virtual arma::mat get_Cygivenx() const=0;
  
  /**
   * @brief Performs the change data operation.
   */
  virtual void change_data()=0;
  /**
   * @brief Performs the change data operation.
   *
   * @param new_data The new data.
   */
  virtual void change_data(std::shared_ptr<Data> new_data)=0;
  
  /**
   * @brief Returns the data.
   *
   * @return The result.
   */
  Data* get_data() const;
  
  /**
   * @brief Returns the measurement pointer.
   *
   * @return The result.
   */
  arma::colvec* get_measurement_pointer();
  /**
   * @brief Returns the measurement pointer.
   *
   * @return The result.
   */
  const arma::colvec* get_measurement_pointer() const;
  
  /**
   * @brief Returns the measurement variables.
   *
   * @return The result.
   */
  std::vector<std::string> get_measurement_variables() const;
  
protected:
  
  friend MeasurementCovarianceEstimatorOutput;
  friend GenericMeasurementCovarianceEstimatorOutput;
  friend DirectGaussianMeasurementCovarianceEstimatorOutput;
  friend GaussianMeasurementCovarianceEstimatorOutput;
  friend MixedGenericDirectGaussianMeasurementCovarianceEstimatorOutput;
  
  // Not stored here. Stored in "main'.
  /** @brief The data. */
  Data* data;
  
  // not stored here
  /** @brief The current data. */
  Data* current_data;
  
  // Not stored here. Stored in "main'.
  /** @brief The rng. */
  RandomNumberGenerator* rng;
  
  // Not stored here. Stored in "main'.
  /** @brief The seed. */
  size_t* seed;
  
  /** @brief The measurement. */
  arma::colvec measurement;
  
  //TransformPtr inverse_transform;
  /** @brief The transform. */
  std::shared_ptr<Transform> transform;
  
  //SummaryStatisticsPtr summary_statistics;
  /** @brief The summary statistics. */
  std::shared_ptr<Transform> summary_statistics;
  
  /**
   * @brief Performs the initialise measurement covariance estimator operation.
   *
   * @return The result.
   */
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator()=0;
  /**
   * @brief Performs the initialise measurement covariance estimator operation.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   *
   * @return The result.
   */
  virtual MeasurementCovarianceEstimatorOutput* initialise_measurement_covariance_estimator(const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Performs the setup measurement variables operation.
   */
  virtual void setup_measurement_variables()=0;
  /**
   * @brief Performs the setup measurement variables operation.
   *
   * @param conditioned_on_parameters The conditioned on parameters.
   */
  virtual void setup_measurement_variables(const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Copies the state of another MeasurementCovarianceEstimator into this object.
   *
   * @param another The MeasurementCovarianceEstimator instance to copy from.
   */
  void make_copy(const MeasurementCovarianceEstimator &another);
  
  /** @brief The set using parameters. */
  bool set_using_parameters;
  
  /** @brief The measurement variables. */
  std::vector<std::string> measurement_variables;
  
};
}

#endif
