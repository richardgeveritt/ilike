#ifndef MATRIXPARAMETERESTIMATOR_H
#define MATRIXPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "parameter_estimator.h"
#include "particles.h"

namespace ilike
{
  /**
   * @file matrix_parameter_estimator.h
   * @brief Defines the MatrixParameterEstimator class.
   *
   * Estimates matrix parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class MatrixParameterEstimator
   * @brief A matrix parameter estimator derived from ParameterEstimator.
   */


class MatrixParameterEstimator : public ParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for MatrixParameterEstimator.
   */
  MatrixParameterEstimator();
  
  /**
   * @brief Destructor for MatrixParameterEstimator.
   */
  virtual ~MatrixParameterEstimator();
  
  /**
   * @brief Copy constructor for MatrixParameterEstimator.
   *
   * @param another The MatrixParameterEstimator instance to copy from.
   */
  MatrixParameterEstimator(const MatrixParameterEstimator &another);
  
  /**
   * @brief Assignment operator for MatrixParameterEstimator.
   *
   * @param another The MatrixParameterEstimator instance to copy from.
   */
  void operator=(const MatrixParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this MatrixParameterEstimator object.
   *
   * @return The result.
   */
  virtual MatrixParameterEstimator* duplicate() const=0;
  
  virtual void fit(const arma::mat &matrix_points,
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
  
  arma::mat estimated;
  
protected:
  
  /**
   * @brief Copies the state of another MatrixParameterEstimator into this object.
   *
   * @param another The MatrixParameterEstimator instance to copy from.
   */
  void make_copy(const MatrixParameterEstimator &another);
  
};
}

#endif
