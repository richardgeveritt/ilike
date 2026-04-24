#ifndef SEQUENTIALDENSITYLIKELIHOODESTIMATORWORKER_H
#define SEQUENTIALDENSITYLIKELIHOODESTIMATORWORKER_H

#include "density_likelihood_estimator_worker.h"

#include <vector>
#include "parameters.h"

namespace ilike
{
  /**
   * @file sequential_density_likelihood_estimator_worker.h
   * @brief Defines the DensityLikelihoodEstimatorOutput class.
   *
   * Stores the output of a density likelihood estimator evaluation. Provides access to log-likelihood values, simulated summaries, and gradient information computed during likelihood estimation.
   *
   * @namespace ilike
   * @class DensityLikelihoodEstimatorOutput
   * @brief The density likelihood estimator output class.
   */


class DensityLikelihoodEstimatorOutput;

class SequentialDensityLikelihoodEstimatorWorker : public DensityLikelihoodEstimatorWorker
{
public:
  
  /**
   * @brief Performs the sequentialdensitylikelihoodestimatorworker operation.
   *
   * @param void The void.
   */
  SequentialDensityLikelihoodEstimatorWorker(void);
  /**
   * @brief Performs the ~sequentialdensitylikelihoodestimatorworker operation.
   *
   * @param void The void.
   */
  virtual ~SequentialDensityLikelihoodEstimatorWorker(void);
  
  /**
   * @brief Performs the sequentialdensitylikelihoodestimatorworker operation.
   *
   * @param the_dle_in The the dle.
   */
  SequentialDensityLikelihoodEstimatorWorker(DensityLikelihoodEstimator* the_dle_in);
  
  /**
   * @brief Performs the sequentialdensitylikelihoodestimatorworker operation.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  SequentialDensityLikelihoodEstimatorWorker(const SequentialDensityLikelihoodEstimatorWorker &another);
  /**
   * @brief Assignment operator for DensityLikelihoodEstimatorOutput.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  void operator=(const SequentialDensityLikelihoodEstimatorWorker &another);
  /**
   * @brief Creates a deep copy of this DensityLikelihoodEstimatorOutput object.
   *
   * @return The result.
   */
  DensityLikelihoodEstimatorWorker* duplicate() const;
  
  /**
   * @brief Returns the points.
   *
   * @return The result.
   */
  std::vector<Parameters> get_points() const;
  
protected:
  
  /**
   * @brief Class-specific implementation for simulate.
   *
   * @param output The output.
   */
  void specific_simulate(DensityLikelihoodEstimatorOutput* output);
  void specific_simulate(DensityLikelihoodEstimatorOutput* output,
                         const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Performs the subsample specific simulate operation.
   *
   * @param output The output.
   */
  void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output);
  void subsample_specific_simulate(DensityLikelihoodEstimatorOutput* output,
                                   const Parameters &conditioned_on_parameters);
  
  /**
   * @brief Copies the state of another DensityLikelihoodEstimatorOutput into this object.
   *
   * @param another The DensityLikelihoodEstimatorOutput instance to copy from.
   */
  void make_copy(const SequentialDensityLikelihoodEstimatorWorker &another);
  
  /** @brief The points. */
  std::vector<Parameters> points;
  
};
}

#endif
