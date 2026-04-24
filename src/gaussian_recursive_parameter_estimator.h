#ifndef GAUSSIANRECURSIVEPARAMETERESTIMATOR_H
#define GAUSSIANRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include <boost/unordered_map.hpp>
#include "recursive_parameter_estimator.h"
#include "gaussian_proposal_info.h"

namespace ilike
{
  /**
   * @file gaussian_recursive_parameter_estimator.h
   * @brief Defines the GaussianRecursiveParameterEstimator class.
   *
   * A recursive parameter estimator for gaussian parameters. Updates estimates online as new observations arrive, for use in adaptive SMC or MCMC algorithms.
   *
   * @namespace ilike
   * @class GaussianRecursiveParameterEstimator
   * @brief A gaussian recursive parameter estimator derived from RecursiveParameterEstimator.
   */


class GaussianRecursiveParameterEstimator : public RecursiveParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for GaussianRecursiveParameterEstimator.
   */
  GaussianRecursiveParameterEstimator();
  
  /**
   * @brief Destructor for GaussianRecursiveParameterEstimator.
   */
  virtual ~GaussianRecursiveParameterEstimator();
  
  /**
   * @brief Copy constructor for GaussianRecursiveParameterEstimator.
   *
   * @param another The GaussianRecursiveParameterEstimator instance to copy from.
   */
  GaussianRecursiveParameterEstimator(const GaussianRecursiveParameterEstimator &another);
  
  /**
   * @brief Assignment operator for GaussianRecursiveParameterEstimator.
   *
   * @param another The GaussianRecursiveParameterEstimator instance to copy from.
   */
  void operator=(const GaussianRecursiveParameterEstimator &another);
  /**
   * @brief Creates a deep copy and returns it as a gaussian pointer.
   *
   * @return The result.
   */
  virtual GaussianRecursiveParameterEstimator* gaussian_duplicate() const=0;
  
  GaussianProposalInfo estimated;
  
protected:
  
  /**
   * @brief Copies the state of another GaussianRecursiveParameterEstimator into this object.
   *
   * @param another The GaussianRecursiveParameterEstimator instance to copy from.
   */
  void make_copy(const GaussianRecursiveParameterEstimator &another);
  
};
}

#endif
