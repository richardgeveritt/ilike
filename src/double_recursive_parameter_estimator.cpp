#include "double_recursive_parameter_estimator.h"

namespace ilike
{
DoubleRecursiveParameterEstimator::DoubleRecursiveParameterEstimator()
:RecursiveParameterEstimator()
{
}

DoubleRecursiveParameterEstimator::DoubleRecursiveParameterEstimator(double initial_value)
{
  this->estimated = initial_value;
}

DoubleRecursiveParameterEstimator::~DoubleRecursiveParameterEstimator()
{
  
}

//Copy constructor for the DoubleRecursiveParameterEstimator class.
DoubleRecursiveParameterEstimator::DoubleRecursiveParameterEstimator(const DoubleRecursiveParameterEstimator &another)
:RecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void DoubleRecursiveParameterEstimator::operator=(const DoubleRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  RecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

void DoubleRecursiveParameterEstimator::make_copy(const DoubleRecursiveParameterEstimator &another)
{
  this->estimated = another.estimated;
}
}
