#include "vector_recursive_parameter_estimator.h"

namespace ilike
{
VectorRecursiveParameterEstimator::VectorRecursiveParameterEstimator()
:RecursiveParameterEstimator()
{
}

VectorRecursiveParameterEstimator::~VectorRecursiveParameterEstimator()
{
  
}

//Copy constructor for the VectorRecursiveParameterEstimator class.
VectorRecursiveParameterEstimator::VectorRecursiveParameterEstimator(const VectorRecursiveParameterEstimator &another)
:RecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void VectorRecursiveParameterEstimator::operator=(const VectorRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  RecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

void VectorRecursiveParameterEstimator::make_copy(const VectorRecursiveParameterEstimator &another)
{
  this->estimated = another.estimated;
}
}
