#include "sample_mean_vector_parameter_estimator.h"
#include "utils.h"

SampleMeanVectorParameterEstimator::SampleMeanVectorParameterEstimator()
  :VectorParameterEstimator()
{
}

SampleMeanVectorParameterEstimator::~SampleMeanVectorParameterEstimator()
{
  
}

//Copy constructor for the SampleMeanVectorParameterEstimator class.
SampleMeanVectorParameterEstimator::SampleMeanVectorParameterEstimator(const SampleMeanVectorParameterEstimator &another)
  :VectorParameterEstimator(another)
{
  this->make_copy(another);
}

void SampleMeanVectorParameterEstimator::operator=(const SampleMeanVectorParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  VectorParameterEstimator::operator=(another);
  this->make_copy(another);
}

VectorParameterEstimator* SampleMeanVectorParameterEstimator::duplicate(void) const
{
  return( new SampleMeanVectorParameterEstimator(*this));
}

void SampleMeanVectorParameterEstimator::make_copy(const SampleMeanVectorParameterEstimator &another)
{
}

void SampleMeanVectorParameterEstimator::fit(const arma::mat &matrix_points,
                                             const arma::colvec &wt)
{
  this->estimated = arma::conv_to< arma::colvec >::from(mean_wt(matrix_points,wt));
}
