#ifndef SMCADAPTOR_H
#define SMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file smc_adaptor.h
   * @brief Defines the SMCOutput class.
   *
   * Stores and manages the output produced by SMC. Holds results such as log-likelihood estimates, samples, or diagnostics returned after running the associated algorithm.
   *
   * @namespace ilike
   * @class SMCOutput
   * @brief The smc output class.
   */


class SMCOutput;
class EnsembleKalmanOutput;

class SMCAdaptor
{
  
public:
  
  /**
   * @brief Performs the smcadaptor operation.
   */
  SMCAdaptor();
  /**
   * @brief Performs the ~smcadaptor operation.
   */
  virtual ~SMCAdaptor();
  
  /**
   * @brief Performs the smcadaptor operation.
   *
   * @param another The SMCOutput instance to copy from.
   */
  SMCAdaptor(const SMCAdaptor &another);
  
  /**
   * @brief Assignment operator for SMCOutput.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void operator=(const SMCAdaptor &another);
  /**
   * @brief Creates a deep copy of this SMCOutput object.
   *
   * @return The result.
   */
  virtual SMCAdaptor* duplicate() const=0;
  
  /**
   * @brief Performs the smc adapt operation.
   *
   * @param current_state The current state.
   */
  virtual void smc_adapt(SMCOutput* current_state)=0;
  /**
   * @brief Performs the ensemble adapt operation.
   *
   * @param current_state The current state.
   */
  virtual void ensemble_adapt(EnsembleKalmanOutput* current_state)=0;
  
  // Mh has its own parameters.
  // Stochastic has some weights.
  // MH has sim prop and eval prop, take in params. Use current value in acceptance, Set current value if accepted.
  // Proposal needs to call simulate in all llhdoutputs
  
protected:
  
  /**
   * @brief Copies the state of another SMCOutput into this object.
   *
   * @param another The SMCOutput instance to copy from.
   */
  void make_copy(const SMCAdaptor &another);
  
};
}

#endif
