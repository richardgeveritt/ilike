#ifndef GAUSSIANMCMCADAPTOR_H
#define GAUSSIANMCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "mcmc_adaptor.h"
#include "particle.h"
#include "distributions.h"
#include "scale.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file gaussian_mcmc_adaptor.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class ScaleRecursiveParameterEstimator;
class VectorRecursiveParameterEstimator;
class GaussianRecursiveParameterEstimator;

class GaussianMCMCAdaptor : public MCMCAdaptor
{
  
public:
  
  /**
   * @brief Performs the gaussianmcmcadaptor operation.
   */
  GaussianMCMCAdaptor();
  /**
   * @brief Performs the ~gaussianmcmcadaptor operation.
   */
  virtual ~GaussianMCMCAdaptor();
  
  /**
   * @brief Performs the gaussianmcmcadaptor operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  GaussianMCMCAdaptor(const GaussianMCMCAdaptor &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const GaussianMCMCAdaptor &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  MCMCAdaptor* duplicate() const;
  
protected:
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const GaussianMCMCAdaptor &another);
  
  void specific_mcmc_adapt(const Particle &latest_particle,
                           size_t iteration_counter);
  
  /** @brief The initial proposal info. */
  GaussianProposalInfo initial_proposal_info;
  
  /** @brief The variable names. */
  std::vector<std::string> variable_names;
  
  // not stored here
  /** @brief The gaussian info pointers. */
  std::vector<GaussianProposalInfo*> gaussian_info_pointers;
  
  // stored here
  /** @brief The scale estimators. */
  std::vector<ScaleRecursiveParameterEstimator*> scale_estimators;
  /** @brief The mean estimators. */
  std::vector<VectorRecursiveParameterEstimator*> mean_estimators;
  /** @brief The gaussian estimators. */
  std::vector<GaussianRecursiveParameterEstimator*> gaussian_estimators;
  
};
}

#endif
