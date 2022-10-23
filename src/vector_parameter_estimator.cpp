#include "vector_parameter_estimator.h"

VectorParameterEstimator::VectorParameterEstimator()
  :ParameterEstimator()
{
}

VectorParameterEstimator::~VectorParameterEstimator()
{
  
}

//Copy constructor for the VectorParameterEstimator class.
VectorParameterEstimator::VectorParameterEstimator(const VectorParameterEstimator &another)
  :ParameterEstimator(another)
{
  this->make_copy(another);
}

void VectorParameterEstimator::operator=(const VectorParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  ParameterEstimator::operator=(another);
  this->make_copy(another);
}

void VectorParameterEstimator::make_copy(const VectorParameterEstimator &another)
{
  this->estimated = another.estimated;
}
