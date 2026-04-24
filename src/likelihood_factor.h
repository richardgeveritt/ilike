#ifndef LIKELIHOODFACTOR_H
#define LIKELIHOODFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "factor.h"

namespace ilike
{
  /**
   * @file likelihood_factor.h
   * @brief Defines the LikelihoodFactor class.
   *
   * A generic likelihood factor. Evaluates the log-likelihood contribution from a generic model for use inside an ExactLikelihoodEstimator.
   *
   * @namespace ilike
   * @class LikelihoodFactor
   * @brief A likelihood factor derived from Factor.
   */


class LikelihoodFactor : public Factor
{
  
public:
  
  /**
   * @brief Default constructor for LikelihoodFactor.
   */
  LikelihoodFactor();
  /**
   * @brief Constructs a LikelihoodFactor object.
   *
   * @param data_in The data.
   */
  LikelihoodFactor(Data* data_in);
  /**
   * @brief Destructor for LikelihoodFactor.
   */
  virtual ~LikelihoodFactor();
  
  /**
   * @brief Copy constructor for LikelihoodFactor.
   *
   * @param another The LikelihoodFactor instance to copy from.
   */
  LikelihoodFactor(const LikelihoodFactor &another);
  
  /**
   * @brief Assignment operator for LikelihoodFactor.
   *
   * @param another The LikelihoodFactor instance to copy from.
   */
  void operator=(const LikelihoodFactor &another);
  /**
   * @brief Creates a deep copy and returns it as a likelihood_factor pointer.
   *
   * @return The result.
   */
  virtual LikelihoodFactor* likelihood_factor_duplicate() const=0;
  
  /**
   * @brief Evaluates.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double evaluate(const Parameters &input) const;
  /**
   * @brief Performs the likelihood evaluate operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual double likelihood_evaluate(const Parameters &input) const=0;
  
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &input) const;
  virtual arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                                 const Parameters &input) const=0;
  
  /**
   * @brief Sets the data.
   *
   * @param data_in The data.
   */
  void set_data(Data* data_in);
  
protected:
  
  /**
   * @brief Class-specific implementation for set data.
   */
  virtual void specific_set_data()=0;
  
  /**
   * @brief Copies the state of another LikelihoodFactor into this object.
   *
   * @param another The LikelihoodFactor instance to copy from.
   */
  void make_copy(const LikelihoodFactor &another);
  
  // Not stored here. Stored in "main'.
  /** @brief The data. */
  Data* data;
  
};
}

#endif
