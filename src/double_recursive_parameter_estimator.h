#ifndef MATRIXPARAMETERESTIMATOR_H
#define MATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "recursive_parameter_estimator.h"

namespace ilike
{
  /**
   * @file double_recursive_parameter_estimator.h
   * @brief Defines the DoubleRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for double parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class DoubleRecursiveParameterEstimator
   * @brief A double recursive parameter estimator derived from RecursiveParameterEstimator.
   */


class DoubleRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for DoubleRecursiveParameterEstimator.
   */
  DoubleRecursiveParameterEstimator();
  /**
   * @brief Constructs a DoubleRecursiveParameterEstimator object.
   *
   * @param initial_value The initial value.
   */
  DoubleRecursiveParameterEstimator(double initial_value);
  
  /**
   * @brief Destructor for DoubleRecursiveParameterEstimator.
   */
  virtual ~DoubleRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for DoubleRecursiveParameterEstimator.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  DoubleRecursiveParameterEstimator(const DoubleRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for DoubleRecursiveParameterEstimator.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const DoubleRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a double pointer.
   *
   * @return The result.
   */
  virtual DoubleRecursiveParameterEstimator* double_duplicate() const=0;
  
  double estimated;
  
protected:
  
  /**
   * @brief Copies the state of another DoubleRecursiveParameterEstimator into this object.
   *
   * @param another The DoubleRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const DoubleRecursiveParameterEstimator &another);
  
};
}

#endif
