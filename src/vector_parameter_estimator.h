#ifndef VECTORPARAMETERESTIMATOR_H
#define VECTORPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

namespace ilike
{
  /**
   * @file vector_parameter_estimator.h
   * @brief Defines the VectorParameterEstimator class.
   *
   * Estimates vector parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class VectorParameterEstimator
   * @brief A vector parameter estimator derived from ParameterEstimator.
   */


class VectorParameterEstimator : public ParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for VectorParameterEstimator.
   */
  VectorParameterEstimator();
  
  /**
   * @brief Destructor for VectorParameterEstimator.
   */
  virtual ~VectorParameterEstimator();
  
  /**
   * @brief Copy constructor for VectorParameterEstimator.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  VectorParameterEstimator(const VectorParameterEstimator &another);
  
  /**
   * @brief Assignment operator for VectorParameterEstimator.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void operator=(const VectorParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this VectorParameterEstimator object.
   *
   * @return The result.
   */
  virtual VectorParameterEstimator* duplicate() const=0;
  
  virtual void fit(const arma::mat &points,
                   const arma::colvec &normalised_log_weights)=0;
  
  void fit(const std::vector<std::string> &variables,
           const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<Parameters> &points,
           const arma::colvec &normalised_log_weights);
  
  void fit(const std::string &variable,
           const std::vector<MoveOutput*> &points,
           const arma::colvec &normalised_log_weights);
  
  arma::colvec estimated;
  
protected:
  
  /**
   * @brief Copies the state of another VectorParameterEstimator into this object.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void make_copy(const VectorParameterEstimator &another);
  
};
}

#endif
