#include "parameter_estimator.h"

ParameterEstimator::ParameterEstimator()
{
}

ParameterEstimator::~ParameterEstimator()
{
}

ParameterEstimator::ParameterEstimator(const ParameterEstimator &another)
{
  this->make_copy(another);
}

void ParameterEstimator::operator=(const ParameterEstimator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void ParameterEstimator::make_copy(const ParameterEstimator &another)
{
}
