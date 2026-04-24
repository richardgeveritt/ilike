#ifndef DENSITYESTIMATEADAPTOR_H
#define DENSITYESTIMATEADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "smc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"

namespace ilike
{
  /**
   * @file density_estimate_adaptor.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;

class DensityEstimateAdaptor : public SMCAdaptor
{
  
public:
  
  /**
   * @brief Performs the densityestimateadaptor operation.
   */
  DensityEstimateAdaptor();
  /**
   * @brief Performs the ~densityestimateadaptor operation.
   */
  virtual ~DensityEstimateAdaptor();
  
  /**
   * @brief Performs the densityestimateadaptor operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  DensityEstimateAdaptor(const DensityEstimateAdaptor &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const DensityEstimateAdaptor &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  SMCAdaptor* duplicate() const;
  
  /**
   * @brief Performs the smc adapt operation.
   *
   * @param current_state The current state.
   */
  void smc_adapt(SMCOutput* current_state);
  /**
   * @brief Performs the ensemble adapt operation.
   *
   * @param current_state The current state.
   */
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  
protected:
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const DensityEstimateAdaptor &another);
  
  // needs density estimator hierarchy
  
};
}

#endif
