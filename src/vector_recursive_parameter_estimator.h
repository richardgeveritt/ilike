#ifndef VECTORRECURSIVEPARAMETERESTIMATOR_H
#define VECTORRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include <boost/unordered_map.hpp>
#include "recursive_parameter_estimator.h"

namespace ilike
{
  /**
   * @file vector_recursive_parameter_estimator.h
   * @brief Defines the VectorRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for vector parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class VectorRecursiveParameterEstimator
   * @brief A vector recursive parameter estimator derived from RecursiveParameterEstimator.
   */


class VectorRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for VectorRecursiveParameterEstimator.
   */
  VectorRecursiveParameterEstimator();
  
  /**
   * @brief Destructor for VectorRecursiveParameterEstimator.
   */
  virtual ~VectorRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for VectorRecursiveParameterEstimator.
   *
   * @param another The VectorRecursiveParameterEstimator instance to copy from.
   */
  VectorRecursiveParameterEstimator(const VectorRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for VectorRecursiveParameterEstimator.
   *
   * @param another The VectorRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const VectorRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a vector pointer.
   *
   * @return The result.
   */
  virtual VectorRecursiveParameterEstimator* vector_duplicate() const=0;
  
  arma::colvec estimated;
  
protected:
  
  /**
   * @brief Copies the state of another VectorRecursiveParameterEstimator into this object.
   *
   * @param another The VectorRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const VectorRecursiveParameterEstimator &another);
  
};
}

#endif
