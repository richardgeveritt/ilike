#ifndef CUSTOMDISTRIBUTIONFACTOR_H
#define CUSTOMDISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file custom_distribution_factor.h
   * @brief Defines the CustomDistributionFactor class.
   *
   * A custom distribution factor. Evaluates the log-density of a custom prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class CustomDistributionFactor
   * @brief A custom distribution factor derived from DistributionFactor.
   */


class CustomDistributionFactor : public DistributionFactor
{
  
public:
  
  /**
   * @brief Default constructor for CustomDistributionFactor.
   */
  CustomDistributionFactor();
  
  /**
   * @brief Constructs a CustomDistributionFactor object.
   *
   * @param distribution_in The distribution.
   */
  CustomDistributionFactor(EvaluateLogDistributionPtr distribution_in);
  
  CustomDistributionFactor(EvaluateLogDistributionPtr distribution_in,
                           EvaluateGradientLogDistributionPtr distribution_gradient_in);
  
  /**
   * @brief Destructor for CustomDistributionFactor.
   */
  virtual ~CustomDistributionFactor();
  
  /**
   * @brief Copy constructor for CustomDistributionFactor.
   *
   * @param another The CustomDistributionFactor instance to copy from.
   */
  CustomDistributionFactor(const CustomDistributionFactor &another);
  
  /**
   * @brief Assignment operator for CustomDistributionFactor.
   *
   * @param another The CustomDistributionFactor instance to copy from.
   */
  void operator=(const CustomDistributionFactor &another);
  /**
   * @brief Creates a deep copy of this CustomDistributionFactor object.
   *
   * @return The result.
   */
  Factor* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a distribution_factor pointer.
   *
   * @return The result.
   */
  DistributionFactor* distribution_factor_duplicate() const;
  
  /**
   * @brief Performs the distribution evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double distribution_evaluate(const Parameters &input) const;
  
  arma::mat distribution_evaluate_gradient(const std::string &variable,
                                           const Parameters &input) const;
  
protected:
  
  /**
   * @brief Copies the state of another CustomDistributionFactor into this object.
   *
   * @param another The CustomDistributionFactor instance to copy from.
   */
  void make_copy(const CustomDistributionFactor &another);
  
  /** @brief The distribution. */
  EvaluateLogDistributionPtr distribution;
  /** @brief The distribution gradient. */
  EvaluateGradientLogDistributionPtr distribution_gradient;
  
};
}

#endif
