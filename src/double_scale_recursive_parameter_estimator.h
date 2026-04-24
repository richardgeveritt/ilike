#ifndef DOUBLESCALERECURSIVEPARAMETERESTIMATOR_H
#define DOUBLESCALERECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "particle.h"
#include "scale_recursive_parameter_estimator.h"

namespace ilike
{
  /**
   * @file double_scale_recursive_parameter_estimator.h
   * @brief Defines the DoubleRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for double parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class DoubleRecursiveParameterEstimator
   * @brief The double recursive parameter estimator class.
   */


class DoubleRecursiveParameterEstimator;

class DoubleScaleRecursiveParameterEstimator : public ScaleRecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Performs the doublescalerecursiveparameterestimator operation.
   */
  DoubleScaleRecursiveParameterEstimator();
  
  /**
   * @brief Performs the ~doublescalerecursiveparameterestimator operation.
   */
  virtual ~DoubleScaleRecursiveParameterEstimator();
  
  /**
   * @brief Performs the doublescalerecursiveparameterestimator operation.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  DoubleScaleRecursiveParameterEstimator(const DoubleScaleRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for DoubleRecursiveParameterEstimator.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const DoubleScaleRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this DoubleRecursiveParameterEstimator object.
   *
   * @return The result.
   */
  RecursiveParameterEstimator* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a scale pointer.
   *
   * @return The result.
   */
  ScaleRecursiveParameterEstimator* scale_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  // stored here
  /** @brief The recursive estimator. */
  DoubleRecursiveParameterEstimator* recursive_estimator;
  
  /**
   * @brief Copies the state of another DoubleRecursiveParameterEstimator into this object.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const DoubleScaleRecursiveParameterEstimator &another);
  
};
}

#endif
