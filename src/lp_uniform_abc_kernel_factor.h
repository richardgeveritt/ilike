#ifndef EUCLIDEANUNIFORMABCKERNELFACTOR_H
#define EUCLIDEANUNIFORMABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "abc_kernel_factor.h"

namespace ilike
{
  /**
   * @file lp_uniform_abc_kernel_factor.h
   * @brief Defines the LpUniformABCKernelFactor class.
   *
   * A lp uniform ABC kernel factor. Evaluates the ABC kernel (a lp uniform function of the distance between observed and simulated data) for use in an ABCLikelihoodEstimator.
   *
   * @namespace ilike
   * @class LpUniformABCKernelFactor
   * @brief A lp uniform abc kernel factor derived from ABCKernelFactor.
   */



class LpUniformABCKernelFactor : public ABCKernelFactor
{
  
public:
  
  /**
   * @brief Default constructor for LpUniformABCKernelFactor.
   */
  LpUniformABCKernelFactor();
  LpUniformABCKernelFactor(double p_in,
                           const std::vector<std::string> &data_variables_in,
                           const std::string &epsilon_variable_in,
                           Data* data_in);
  LpUniformABCKernelFactor(double p_in,
                           const std::vector<std::string> &data_variables_in,
                           const std::string &epsilon_variable_in,
                           const std::string &scale_variable_in,
                           Data* data_in);
  
  /**
   * @brief Destructor for LpUniformABCKernelFactor.
   */
  virtual ~LpUniformABCKernelFactor();
  
  /**
   * @brief Copy constructor for LpUniformABCKernelFactor.
   *
   * @param another The LpUniformABCKernelFactor instance to copy from.
   */
  LpUniformABCKernelFactor(const LpUniformABCKernelFactor &another);
  
  /**
   * @brief Assignment operator for LpUniformABCKernelFactor.
   *
   * @param another The LpUniformABCKernelFactor instance to copy from.
   */
  void operator=(const LpUniformABCKernelFactor &another);
  /**
   * @brief Creates a deep copy of this LpUniformABCKernelFactor object.
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
   * @brief Copies the state of another LpUniformABCKernelFactor into this object.
   *
   * @param another The LpUniformABCKernelFactor instance to copy from.
   */
  void make_copy(const LpUniformABCKernelFactor &another);
  
  /** @brief The constant. */
  double constant;
  
  /** @brief The p. */
  double p;
  
};

}

#endif
