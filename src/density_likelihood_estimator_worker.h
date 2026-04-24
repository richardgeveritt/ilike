#ifndef DENSITYLIKELIHOODESTIMATORWORKER_H
#define DENSITYLIKELIHOODESTIMATORWORKER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "distributions.h"
#include "parameters.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file density_likelihood_estimator_worker.h
   * @brief Defines the DensityLikelihoodEstimator class.
   *
   * Estimates the density likelihood for a given set of parameters. Implements the estimator interface for use inside SMC, MCMC, or importance-sampling algorithms.
   *
   * @namespace ilike
   * @class DensityLikelihoodEstimator
   * @brief The density likelihood estimator class.
   */


class DensityLikelihoodEstimator;
class DensityLikelihoodEstimatorOutput;

class DensityLikelihoodEstimatorWorker
{
public:
  
  /**
   * @brief Performs the densitylikelihoodestimatorworker operation.
   */
  DensityLikelihoodEstimatorWorker();
  
  /**
   * @brief Performs the densitylikelihoodestimatorworker operation.
   *
   * @param the_dle_in The the dle.
   */
  DensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in);
  
  /**
   * @brief Performs the ~densitylikelihoodestimatorworker operation.
   *
   * @param void The void.
   */
  virtual ~DensityLikelihoodEstimatorWorker(void);
  
  /**
   * @brief Performs the densitylikelihoodestimatorworker operation.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  DensityLikelihoodEstimatorWorker(const DensityLikelihoodEstimatorWorker &another);
  /**
   * @brief Assignment operator for DensityLikelihoodEstimator.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  void operator=(const DensityLikelihoodEstimatorWorker &another);
  /**
   * @brief Creates a deep copy of this DensityLikelihoodEstimator object.
   *
   * @return The result.
   */
  virtual DensityLikelihoodEstimatorWorker* duplicate() const=0;
  
  /**
   * @brief Returns the points.
   *
   * @return The result.
   */
  virtual std::vector<Parameters> get_points() const=0;
  /**
   * @brief Returns the number of points.
   *
   * @return The result.
   */
  size_t get_number_of_points() const;
  /**
   * @brief Returns the rng.
   *
   * @return The result.
   */
  RandomNumberGenerator* get_rng();
  /**
   * @brief Returns the seed.
   *
   * @return The result.
   */
  size_t get_seed() const;
  /**
   * @brief Sets the seed.
   *
   * @param seed_in The seed.
   */
  void set_seed(size_t seed_in);
  
  /**
   * @brief Simulates the required variables.
   *
   * @param output The output.
   */
  void simulate(DensityLikelihoodEstimatorOutput* output);
  void simulate(DensityLikelihoodEstimatorOutput* output,
                const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the subsample simulate operation.
   *
   * @param output The output.
   */
  void subsample_simulate(DensityLikelihoodEstimatorOutput* output);
  void subsample_simulate(DensityLikelihoodEstimatorOutput* output,
                          const Parameters &conditioned_on_parameters);
  
protected:
  
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param output The output.
   */
  virtual void specific_simulate(DensityLikelihoodEstimatorOutput* output)=0;
  virtual void specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                 const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Performs the subsample specific simulate operation.
   *
   * @param output The output.
   */
  virtual void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output)=0;
  virtual void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                           const Parameters &conditioned_on_parameters)=0;
  
  /**
   * @brief Copies the state of another DensityLikelihoodEstimator into this object.
   *
   * @param another The DensityLikelihoodEstimator instance to copy from.
   */
  void make_copy(const DensityLikelihoodEstimatorWorker &another);
  
  /** @brief The the dle. */
  DensityLikelihoodEstimator* the_dle; // not stored here
  
};
}

#endif
