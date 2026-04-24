#ifndef PARAMETERESTIMATOR_H
#define PARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

namespace ilike
{
  /**
   * @file parameter_estimator.h
   * @brief Defines the ParameterEstimator class.
   *
   * Estimates generic parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class ParameterEstimator
   * @brief The parameter estimator class.
   */


class ParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for ParameterEstimator.
   */
  ParameterEstimator();
  /**
   * @brief Destructor for ParameterEstimator.
   */
  virtual ~ParameterEstimator();
  
  /**
   * @brief Copy constructor for ParameterEstimator.
   *
   * @param another The ParameterEstimator instance to copy from.
   */
  ParameterEstimator(const ParameterEstimator &another);
  
  /**
   * @brief Assignment operator for ParameterEstimator.
   *
   * @param another The ParameterEstimator instance to copy from.
   */
  void operator=(const ParameterEstimator &another);
  
protected:
  
  /**
   * @brief Copies the state of another ParameterEstimator into this object.
   *
   * @param another The ParameterEstimator instance to copy from.
   */
  void make_copy(const ParameterEstimator &another);
  
};
}

#endif
