#ifndef SAMPLEMEANVECTORPARAMETERESTIMATOR_H
#define SAMPLEMEANVECTORPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "vector_parameter_estimator.h"

namespace ilike
{
  /**
   * @file sample_mean_vector_parameter_estimator.h
   * @brief Defines the SampleMeanVectorParameterEstimator class.
   *
   * Estimates sample mean vector parameters from a particle set or ensemble. Used in adaptive algorithms to tune proposal or model parameters.
   *
   * @namespace ilike
   * @class SampleMeanVectorParameterEstimator
   * @brief A sample mean vector parameter estimator derived from VectorParameterEstimator.
   */


class SampleMeanVectorParameterEstimator : public VectorParameterEstimator
{
  
public:
  
  /**
   * @brief Default constructor for SampleMeanVectorParameterEstimator.
   */
  SampleMeanVectorParameterEstimator();
  
  /**
   * @brief Destructor for SampleMeanVectorParameterEstimator.
   */
  virtual ~SampleMeanVectorParameterEstimator();
  
  /**
   * @brief Copy constructor for SampleMeanVectorParameterEstimator.
   *
   * @param another The SampleMeanVectorParameterEstimator instance to copy from.
   */
  SampleMeanVectorParameterEstimator(const SampleMeanVectorParameterEstimator &another);
  
  /**
   * @brief Assignment operator for SampleMeanVectorParameterEstimator.
   *
   * @param another The SampleMeanVectorParameterEstimator instance to copy from.
   */
  void operator=(const SampleMeanVectorParameterEstimator &another);
  /**
   * @brief Creates a deep copy of this SampleMeanVectorParameterEstimator object.
   *
   * @return The result.
   */
  VectorParameterEstimator* duplicate() const;
  
  void fit(const arma::mat &points,
           const arma::colvec &normalised_log_weights);
  
protected:
  
  /**
   * @brief Copies the state of another SampleMeanVectorParameterEstimator into this object.
   *
   * @param another The SampleMeanVectorParameterEstimator instance to copy from.
   */
  void make_copy(const SampleMeanVectorParameterEstimator &another);
  
};
}

#endif
