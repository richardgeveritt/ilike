#include "scale_recursive_parameter_estimator.h"

ScaleRecursiveParameterEstimator::ScaleRecursiveParameterEstimator()
  :RecursiveParameterEstimator()
{
}

ScaleRecursiveParameterEstimator::~ScaleRecursiveParameterEstimator()
{
  
}

//Copy constructor for the ScaleRecursiveParameterEstimator class.
ScaleRecursiveParameterEstimator::ScaleRecursiveParameterEstimator(const ScaleRecursiveParameterEstimator &another)
  :RecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void ScaleRecursiveParameterEstimator::operator=(const ScaleRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  RecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

void ScaleRecursiveParameterEstimator::make_copy(const ScaleRecursiveParameterEstimator &another)
{
  this->estimated = another.estimated;
}
