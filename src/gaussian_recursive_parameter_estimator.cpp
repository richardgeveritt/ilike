#include "gaussian_recursive_parameter_estimator.h"

namespace ilike
{
GaussianRecursiveParameterEstimator::GaussianRecursiveParameterEstimator()
:RecursiveParameterEstimator()
{
}

GaussianRecursiveParameterEstimator::~GaussianRecursiveParameterEstimator()
{
  
}

//Copy constructor for the GaussianRecursiveParameterEstimator class.
GaussianRecursiveParameterEstimator::GaussianRecursiveParameterEstimator(const GaussianRecursiveParameterEstimator &another)
:RecursiveParameterEstimator(another)
{
  this->make_copy(another);
}

void GaussianRecursiveParameterEstimator::operator=(const GaussianRecursiveParameterEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  RecursiveParameterEstimator::operator=(another);
  this->make_copy(another);
}

void GaussianRecursiveParameterEstimator::make_copy(const GaussianRecursiveParameterEstimator &another)
{
  this->estimated = another.estimated;
}
}
