#ifndef SAMPLEMEANVECTORPARAMETERESTIMATOR_H
#define SAMPLEMEANVECTORPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "vector_parameter_estimator.h"

namespace ilike
{
class SampleMeanVectorParameterEstimator : public VectorParameterEstimator
{
  
public:
  
  SampleMeanVectorParameterEstimator();
  
  virtual ~SampleMeanVectorParameterEstimator();
  
  SampleMeanVectorParameterEstimator(const SampleMeanVectorParameterEstimator &another);
  
  void operator=(const SampleMeanVectorParameterEstimator &another);
  VectorParameterEstimator* duplicate() const;
  
  void fit(const arma::mat &points,
           const arma::colvec &normalised_log_weights);
  
protected:
  
  void make_copy(const SampleMeanVectorParameterEstimator &another);
  
};
}

#endif
