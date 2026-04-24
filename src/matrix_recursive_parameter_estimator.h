#ifndef MATRIXRECURSIVEPARAMETERESTIMATOR_H
#define MATRIXRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "recursive_parameter_estimator.h"

namespace ilike
{
  /**
   * @file matrix_recursive_parameter_estimator.h
   * @brief Defines the MatrixRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for matrix parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class MatrixRecursiveParameterEstimator
   * @brief A matrix recursive parameter estimator derived from RecursiveParameterEstimator.
   */


class MatrixRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for MatrixRecursiveParameterEstimator.
   */
  MatrixRecursiveParameterEstimator();
  
  /**
   * @brief Destructor for MatrixRecursiveParameterEstimator.
   */
  virtual ~MatrixRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for MatrixRecursiveParameterEstimator.
   *
   * @param another The MatrixRecursiveParameterEstimator instance to copy from.
   */
  MatrixRecursiveParameterEstimator(const MatrixRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for MatrixRecursiveParameterEstimator.
   *
   * @param another The MatrixRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const MatrixRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a matrix pointer.
   *
   * @return The result.
   */
  virtual MatrixRecursiveParameterEstimator* matrix_duplicate() const=0;
  
  arma::mat estimated;
  
protected:
  
  /**
   * @brief Copies the state of another MatrixRecursiveParameterEstimator into this object.
   *
   * @param another The MatrixRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const MatrixRecursiveParameterEstimator &another);
  
};
}

#endif
