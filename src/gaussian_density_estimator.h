#ifndef GAUSSIANDENSITYESTIMATOR_H
#define GAUSSIANDENSITYESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "density_estimator.h"

namespace ilike
{
  /**
   * @file gaussian_density_estimator.h
   * @brief Defines the VectorParameterEstimator class.
   *
   * Estimates vector parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class VectorParameterEstimator
   * @brief The vector parameter estimator class.
   */


class VectorParameterEstimator;
class MatrixParameterEstimator;
class DensityEstimatorOutput;

class GaussianDensityEstimator : public DensityEstimator
{

public:

  /**
   * @brief Performs the gaussiandensityestimator operation.
   */
  GaussianDensityEstimator();
  
  /**
   * @brief Performs the gaussiandensityestimator operation.
   *
   * @param variables_in The variables.
   */
  GaussianDensityEstimator(const std::vector<std::string> &variables_in);
  
  GaussianDensityEstimator(const std::vector<std::string> &variables_in,
                           bool unbiased_in);

  /**
   * @brief Performs the ~gaussiandensityestimator operation.
   */
  virtual ~GaussianDensityEstimator();

  /**
   * @brief Performs the gaussiandensityestimator operation.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  GaussianDensityEstimator(const GaussianDensityEstimator &another);

  /**
   * @brief Assignment operator for VectorParameterEstimator.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void operator=(const GaussianDensityEstimator &another);
  /**
   * @brief Creates a deep copy of this VectorParameterEstimator object.
   *
   * @return The result.
   */
  DensityEstimator* duplicate() const;
  
  /**
   * @brief Initialises the estimator and returns an output object.
   *
   * @return The result.
   */
  DensityEstimatorOutput* initialise();
  
  /**
   * @brief Returns the unbiased.
   *
   * @return The result.
   */
  bool get_unbiased() const;

  /*
  void fit(const std::vector<Parameters> &points,
           arma::colvec normalised_log_weights);
  
  double evaluate(const Data &point) const;
  */
  
protected:
  
  // Stored here.
  //VectorParameterEstimator* mean_estimator;
  //MatrixParameterEstimator* covariance_estimator;
  
  /** @brief The unbiased. */
  bool unbiased;

  /**
   * @brief Copies the state of another VectorParameterEstimator into this object.
   *
   * @param another The VectorParameterEstimator instance to copy from.
   */
  void make_copy(const GaussianDensityEstimator &another);

};
}

#endif
