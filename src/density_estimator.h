#ifndef DENSITYESTIMATOR_H
#define DENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

namespace ilike
{
  /**
   * @file density_estimator.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


//class SMCOutput;
class DensityEstimatorOutput;

class DensityEstimator
{
  
public:
  
  /**
   * @brief Performs the densityestimator operation.
   */
  DensityEstimator();
  /**
   * @brief Performs the densityestimator operation.
   *
   * @param variables_in The variables.
   */
  DensityEstimator(const std::vector<std::string> &variables_in);
  /**
   * @brief Performs the ~densityestimator operation.
   */
  virtual ~DensityEstimator();
  
  /**
   * @brief Performs the densityestimator operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  DensityEstimator(const DensityEstimator &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const DensityEstimator &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  virtual DensityEstimator* duplicate() const=0;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  virtual DensityEstimatorOutput* initialise()=0;
  
  /**
   * @brief Returns the variables.
   *
   * @return The result.
   */
  std::vector<std::string> get_variables() const;
  
protected:
  
  /** @brief The variables. */
  std::vector<std::string> variables;
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const DensityEstimator &another);
  
};
}

#endif
