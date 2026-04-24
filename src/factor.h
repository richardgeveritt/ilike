#ifndef FACTOR_H
#define FACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{
  /**
   * @file factor.h
   * @brief Defines the Factor class.
   *
   * Provides factor functionality.
   *
   * @namespace ilike
   * @class Factor
   * @brief The factor class.
   */


class Factor
{
  
public:
  
  /**
   * @brief Default constructor for Factor.
   */
  Factor();
  /**
   * @brief Destructor for Factor.
   */
  virtual ~Factor();
  
  /**
   * @brief Copy constructor for Factor.
   *
   * @param another The Factor instance to copy from.
   */
  Factor(const Factor &another);
  
  /**
   * @brief Assignment operator for Factor.
   *
   * @param another The Factor instance to copy from.
   */
  void operator=(const Factor &another);
  /**
   * @brief Creates a deep copy of this Factor object.
   *
   * @return The result.
   */
  virtual Factor* duplicate() const=0;
  
  /**
   * @brief Evaluates.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual double evaluate(const Parameters &input) const=0;
  virtual arma::mat evaluate_gradient(const std::string &variable,
                                      const Parameters &input) const=0;
  
protected:
  
  /**
   * @brief Copies the state of another Factor into this object.
   *
   * @param another The Factor instance to copy from.
   */
  void make_copy(const Factor &another);
  
};
}

#endif
