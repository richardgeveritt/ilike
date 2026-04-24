#ifndef CUSTOMLIKELIHOODFACTOR_H
#define CUSTOMLIKELIHOODFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ilike_header.h"
#include "likelihood_factor.h"

namespace ilike
{
  /**
   * @file custom_likelihood_factor.h
   * @brief Defines the CustomLikelihoodFactor class.
   *
   * A custom likelihood factor. Evaluates the log-likelihood contribution from a custom model for use inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class CustomLikelihoodFactor
   * @brief A custom likelihood factor derived from LikelihoodFactor.
   */


class CustomLikelihoodFactor : public LikelihoodFactor
{
  
public:
  
  /**
   * @brief Default constructor for CustomLikelihoodFactor.
   */
  CustomLikelihoodFactor();
  CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                         Data* data_in);
  CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                         EvaluateGradientLogLikelihoodPtr likelihood_gradient_in,
                         Data* data_in);
  
  /**
   * @brief Destructor for CustomLikelihoodFactor.
   */
  virtual ~CustomLikelihoodFactor();
  
  /**
   * @brief Copy constructor for CustomLikelihoodFactor.
   *
   * @param another The CustomLikelihoodFactor instance to copy from.
   */
  CustomLikelihoodFactor(const CustomLikelihoodFactor &another);
  
  /**
   * @brief Assignment operator for CustomLikelihoodFactor.
   *
   * @param another The CustomLikelihoodFactor instance to copy from.
   */
  void operator=(const CustomLikelihoodFactor &another);
  /**
   * @brief Creates a deep copy of this CustomLikelihoodFactor object.
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
   * @brief Performs the likelihood evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double likelihood_evaluate(const Parameters &input) const;
  
  arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                         const Parameters &input) const;
  
protected:
  
  /**
   * @brief Class-specific implementation for set data.
   */
  void specific_set_data();
  
  /**
   * @brief Copies the state of another CustomLikelihoodFactor into this object.
   *
   * @param another The CustomLikelihoodFactor instance to copy from.
   */
  void make_copy(const CustomLikelihoodFactor &another);
  
  /** @brief The likelihood. */
  EvaluateLogLikelihoodPtr likelihood;
  /** @brief The likelihood gradient. */
  EvaluateGradientLogLikelihoodPtr likelihood_gradient;
  
};
}

#endif
