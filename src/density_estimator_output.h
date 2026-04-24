#ifndef DENSITYESTIMATOROUTPUT_H
#define DENSITYESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "parameters.h"

//class SMCOutput;

namespace ilike
{
  /**
   * @file density_estimator_output.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class DensityEstimatorOutput
{
  
public:
  
  /**
   * @brief Performs the densityestimatoroutput operation.
   */
  DensityEstimatorOutput();
  /**
   * @brief Performs the ~densityestimatoroutput operation.
   */
  virtual ~DensityEstimatorOutput();
  
  /**
   * @brief Performs the densityestimatoroutput operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  DensityEstimatorOutput(const DensityEstimatorOutput &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const DensityEstimatorOutput &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  virtual DensityEstimatorOutput* duplicate() const=0;
  
  virtual void fit(const std::vector<Parameters> &points,
                   const arma::colvec &normalised_log_weights)=0;
  /**
   * @brief Performs the fit operation.
   *
   * @param points The points.
   */
  void fit(const std::vector<Parameters> &points);
  
  /**
   * @brief Evaluates.
   *
   * @param point The point.
   *
   * @return The result.
   */
  virtual double evaluate(const Parameters &point) const=0;
  
protected:
  
  /** @brief The n. */
  size_t n;
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const DensityEstimatorOutput &another);
  
};
}

#endif
