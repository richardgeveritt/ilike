#ifndef DISTRIBUTIONFACTOR_H
#define DISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "factor.h"

namespace ilike
{
  /**
   * @file distribution_factor.h
   * @brief Defines the DistributionFactor class.
   *
   * A distribution distribution factor. Evaluates the log-density of a distribution prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class DistributionFactor
   * @brief A distribution factor derived from Factor.
   */


class DistributionFactor : public Factor
{
  
public:
  
  /**
   * @brief Default constructor for DistributionFactor.
   */
  DistributionFactor();
  /**
   * @brief Destructor for DistributionFactor.
   */
  virtual ~DistributionFactor();
  
  /**
   * @brief Copy constructor for DistributionFactor.
   *
   * @param another The DistributionFactor instance to copy from.
   */
  DistributionFactor(const DistributionFactor &another);
  
  /**
   * @brief Assignment operator for DistributionFactor.
   *
   * @param another The DistributionFactor instance to copy from.
   */
  void operator=(const DistributionFactor &another);
  /**
   * @brief Creates a deep copy and returns it as a distribution_factor pointer.
   *
   * @return The result.
   */
  virtual DistributionFactor* distribution_factor_duplicate() const=0;
  
  /**
   * @brief Evaluates.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double evaluate(const Parameters &input) const;
  /**
   * @brief Performs the distribution evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual double distribution_evaluate(const Parameters &input) const=0;
  
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &input) const;
  virtual arma::mat distribution_evaluate_gradient(const std::string &variable,
                                                   const Parameters &input) const=0;
  
protected:
  
  /**
   * @brief Copies the state of another DistributionFactor into this object.
   *
   * @param another The DistributionFactor instance to copy from.
   */
  void make_copy(const DistributionFactor &another);
  
};
}

#endif
