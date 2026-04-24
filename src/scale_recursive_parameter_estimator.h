#ifndef SCALERECURSIVEPARAMETERESTIMATOR_H
#define SCALERECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "scale.h"
#include "recursive_parameter_estimator.h"

namespace ilike
{
  /**
   * @file scale_recursive_parameter_estimator.h
   * @brief Defines the ScaleRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for scale parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class ScaleRecursiveParameterEstimator
   * @brief A scale recursive parameter estimator derived from RecursiveParameterEstimator.
   */


class ScaleRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for ScaleRecursiveParameterEstimator.
   */
  ScaleRecursiveParameterEstimator();
  
  /**
   * @brief Destructor for ScaleRecursiveParameterEstimator.
   */
  virtual ~ScaleRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for ScaleRecursiveParameterEstimator.
   *
   * @param another The ScaleRecursiveParameterEstimator instance to copy from.
   */
  ScaleRecursiveParameterEstimator(const ScaleRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for ScaleRecursiveParameterEstimator.
   *
   * @param another The ScaleRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const ScaleRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a scale pointer.
   *
   * @return The result.
   */
  virtual ScaleRecursiveParameterEstimator* scale_duplicate() const=0;
  
  Scale estimated;
  
protected:
  
  /**
   * @brief Copies the state of another ScaleRecursiveParameterEstimator into this object.
   *
   * @param another The ScaleRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const ScaleRecursiveParameterEstimator &another);
  
};
}

#endif
