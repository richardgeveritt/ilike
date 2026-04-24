#ifndef RECURSIVEPARAMETERESTIMATOR_H
#define RECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particle.h"

namespace ilike
{
  /**
   * @file recursive_parameter_estimator.h
   * @brief Defines the ProposalKernel class.
   *
   * A generic proposal kernel. Proposes new parameter values during MCMC or SMC moves using a generic distribution centred on the current state.
   *
   * @namespace ilike
   * @class ProposalKernel
   * @brief The proposal kernel class.
   */


class ProposalKernel;

class RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Performs the recursiveparameterestimator operation.
   */
  RecursiveParameterEstimator();
  /**
   * @brief Performs the ~recursiveparameterestimator operation.
   */
  virtual ~RecursiveParameterEstimator();
  
  /**
   * @brief Performs the recursiveparameterestimator operation.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  RecursiveParameterEstimator(const RecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for ProposalKernel.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void operator=(const RecursiveParameterEstimator &another);
  
  virtual void update(const std::string &variable_name,
                      const Particle &latest_particle,
                      size_t iteration_counter,
                      ProposalKernel* proposal)=0;
  
protected:
  
  /**
   * @brief Copies the state of another ProposalKernel into this object.
   *
   * @param another The ProposalKernel instance to copy from.
   */
  void make_copy(const RecursiveParameterEstimator &another);
  
};

}

#endif
