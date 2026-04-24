#ifndef SCALE_H
#define SCALE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file scale.h
   * @brief Defines the Scale class.
   *
   * Provides scale functionality.
   *
   * @namespace ilike
   * @class Scale
   * @brief The scale class.
   */


class Scale
{
public:
  
  /**
   * @brief Default constructor for Scale.
   */
  Scale();
  /**
   * @brief Destructor for Scale.
   */
  ~Scale();
  
  /**
   * @brief Constructs a Scale object.
   *
   * @param constant_in The constant.
   */
  Scale(double constant_in);
  
  Scale(double constant_in,
        size_t dimension_in);
  
  /**
   * @brief Copy constructor for Scale.
   *
   * @param another The Scale instance to copy from.
   */
  Scale(const Scale &another);
  /**
   * @brief Assignment operator for Scale.
   *
   * @param another The Scale instance to copy from.
   */
  Scale& operator=(const Scale &another);
  
  /**
   * @brief Copy constructor for Scale.
   *
   * @param another The Scale instance to copy from.
   */
  Scale(Scale &&another);
  /**
   * @brief Assignment operator for Scale.
   *
   * @param another The Scale instance to copy from.
   */
  Scale& operator=(Scale &&another);
  
  /**
   * @brief Performs the operator() operation.
   *
   * @return The result.
   */
  double operator()() const;
  
  /**
   * @brief Returns the constant.
   *
   * @return The result.
   */
  double get_constant() const;
  /**
   * @brief Returns the constant.
   *
   * @return The result.
   */
  double& get_constant();
  
  //double operator()(size_t dimension);
  
private:
  
  /**
   * @brief Copies the state of another Scale into this object.
   *
   * @param another The Scale instance to copy from.
   */
  void make_copy(const Scale &another);
  /**
   * @brief Copies the state of another Scale into this object.
   *
   * @param another The Scale instance to copy from.
   */
  void make_copy(Scale &&another);
  
  /** @brief The constant. */
  double constant;
  /** @brief The divide by dimension. */
  bool divide_by_dimension;
  /** @brief The dimension. */
  double dimension;
};
}

#endif
