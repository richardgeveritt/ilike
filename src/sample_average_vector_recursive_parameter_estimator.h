#ifndef SAMPLEAVERAGEVECTORRECURSIVEPARAMETERESTIMATOR_H
#define SAMPLEAVERAGEVECTORRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "vector_recursive_parameter_estimator.h"
#include "ilike_header.h"
#include "particle.h"

namespace ilike
{
  /**
   * @file sample_average_vector_recursive_parameter_estimator.h
   * @brief Defines the SampleAverageVectorRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for sample average vector parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class SampleAverageVectorRecursiveParameterEstimator
   * @brief A sample average vector recursive parameter estimator derived from VectorRecursiveParameterEstimator.
   */


class SampleAverageVectorRecursiveParameterEstimator : public VectorRecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for SampleAverageVectorRecursiveParameterEstimator.
   */
  SampleAverageVectorRecursiveParameterEstimator();
  
  /**
   * @brief Destructor for SampleAverageVectorRecursiveParameterEstimator.
   */
  virtual ~SampleAverageVectorRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for SampleAverageVectorRecursiveParameterEstimator.
   *
   * @param another The SampleAverageVectorRecursiveParameterEstimator instance to copy from.
   */
  SampleAverageVectorRecursiveParameterEstimator(const SampleAverageVectorRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for SampleAverageVectorRecursiveParameterEstimator.
   *
   * @param another The SampleAverageVectorRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const SampleAverageVectorRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this SampleAverageVectorRecursiveParameterEstimator object.
   *
   * @return The result.
   */
  RecursiveParameterEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a vector pointer.
   *
   * @return The result.
   */
  VectorRecursiveParameterEstimator* vector_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  /** @brief The gain. */
  GainPtr gain;
  
  /**
   * @brief Copies the state of another SampleAverageVectorRecursiveParameterEstimator into this object.
   *
   * @param another The SampleAverageVectorRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const SampleAverageVectorRecursiveParameterEstimator &another);
  
};
}

#endif
