#ifndef GAMMADISTRIBUTIONFACTOR_H
#define GAMMADISTRIBUTIONFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distribution_factor.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file gamma_distribution_factor.h
   * @brief Defines the GammaDistributionFactor class.
   *
   * A gamma distribution factor. Evaluates the log-density of a gamma prior as part of the overall model likelihood inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class GammaDistributionFactor
   * @brief A gamma distribution factor derived from DistributionFactor.
   */


class GammaDistributionFactor : public DistributionFactor
{
  
public:
  
  /**
   * @brief Default constructor for GammaDistributionFactor.
   */
  GammaDistributionFactor();
  
  GammaDistributionFactor(const std::string &variable_in,
                          double shape_in,
                          double rate_in);
  
  /**
   * @brief Destructor for GammaDistributionFactor.
   */
  virtual ~GammaDistributionFactor();
  
  /**
   * @brief Copy constructor for GammaDistributionFactor.
   *
   * @param another The GammaDistributionFactor instance to copy from.
   */
  GammaDistributionFactor(const GammaDistributionFactor &another);
  
  /**
   * @brief Assignment operator for GammaDistributionFactor.
   *
   * @param another The GammaDistributionFactor instance to copy from.
   */
  void operator=(const GammaDistributionFactor &another);
  /**
   * @brief Creates a deep copy of this GammaDistributionFactor object.
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
   * @brief Copies the state of another GammaDistributionFactor into this object.
   *
   * @param another The GammaDistributionFactor instance to copy from.
   */
  void make_copy(const GammaDistributionFactor &another);
  
  /** @brief The variable. */
  std::string variable;
  /** @brief The shape. */
  double shape;
  /** @brief The rate. */
  double rate;
  
};
}

#endif
