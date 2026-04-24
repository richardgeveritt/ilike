#ifndef GAUSSIANSMCADAPTOR_H
#define GAUSSIANSMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "smc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file gaussian_smc_adaptor.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class VectorParameterEstimator;
class MatrixParameterEstimator;

class GaussianSMCAdaptor : public SMCAdaptor
{
  
public:
  
  /**
   * @brief Performs the gaussiansmcadaptor operation.
   */
  GaussianSMCAdaptor();
  /**
   * @brief Performs the ~gaussiansmcadaptor operation.
   */
  virtual ~GaussianSMCAdaptor();
  
  /**
   * @brief Performs the gaussiansmcadaptor operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  GaussianSMCAdaptor(const GaussianSMCAdaptor &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const GaussianSMCAdaptor &another);
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
  void make_copy(const GaussianSMCAdaptor &another);
  
  /** @brief The variable names. */
  std::vector<std::string> variable_names;
  
  // not stored here
  /** @brief The gaussian info pointers. */
  std::vector<GaussianProposalInfo*> gaussian_info_pointers;
  
  // stored here
  /** @brief The scales. */
  std::vector<Scale> scales;
  /** @brief The mean estimators. */
  std::vector<VectorParameterEstimator*> mean_estimators;
  /** @brief The covariance estimators. */
  std::vector<MatrixParameterEstimator*> covariance_estimators;
  
};
}

#endif
