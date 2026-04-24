#ifndef GAUSSIANABCKERNELFACTOR_H
#define GAUSSIANABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "abc_kernel_factor.h"

namespace ilike
{
  /**
   * @file gaussian_abc_kernel_factor.h
   * @brief Defines the GaussianABCKernelFactor class.
   *
   * A gaussian ABC kernel factor. Evaluates the ABC kernel (a gaussian function of the distance between observed and simulated data) for use in an ABCLikelihoodEstimator.
   *
   * @namespace ilike
   * @class GaussianABCKernelFactor
   * @brief A gaussian abc kernel factor derived from ABCKernelFactor.
   */



class GaussianABCKernelFactor : public ABCKernelFactor
{
  
public:
  
  /**
   * @brief Default constructor for GaussianABCKernelFactor.
   */
  GaussianABCKernelFactor();
  GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                          const std::string &epsilon_variable_in,
                          Data* data_in);
  GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                          const std::string &epsilon_variable_in,
                          const std::string &weight_variable_in,
                          Data* data_in);
  
  /**
   * @brief Destructor for GaussianABCKernelFactor.
   */
  virtual ~GaussianABCKernelFactor();
  
  /**
   * @brief Copy constructor for GaussianABCKernelFactor.
   *
   * @param another The GaussianABCKernelFactor instance to copy from.
   */
  GaussianABCKernelFactor(const GaussianABCKernelFactor &another);
  
  /**
   * @brief Assignment operator for GaussianABCKernelFactor.
   *
   * @param another The GaussianABCKernelFactor instance to copy from.
   */
  void operator=(const GaussianABCKernelFactor &another);
  /**
   * @brief Creates a deep copy of this GaussianABCKernelFactor object.
   *
   * @return The result.
   */
  Factor* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a likelihood_factor pointer.
   *
   * @return The result.
   */
  LikelihoodFactor* likelihood_factor_duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a abc_kernel_factor pointer.
   *
   * @return The result.
   */
  ABCKernelFactor* abc_kernel_factor_duplicate() const;
  
  /**
   * @brief Performs the likelihood evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double likelihood_evaluate(const Parameters &input) const;
  
  void find_distance(const Parameters &input,
                     double &distance,
                     double &scale_constant) const;
  
  double evaluate_kernel_given_distance(const Parameters &input,
                                        double distance,
                                        double scale_constant) const;
  
  arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                         const Parameters &input) const;
  
protected:
  
  /**
   * @brief Class-specific implementation for set data.
   */
  void specific_set_data();
  
  /**
   * @brief Copies the state of another GaussianABCKernelFactor into this object.
   *
   * @param another The GaussianABCKernelFactor instance to copy from.
   */
  void make_copy(const GaussianABCKernelFactor &another);
  
  /** @brief The constant. */
  double constant;
  
};
}

#endif
