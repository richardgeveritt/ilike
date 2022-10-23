#include "recursive_parameter_estimator.h"
#include "proposal_kernel.h"

RecursiveParameterEstimator::RecursiveParameterEstimator()
{
}

RecursiveParameterEstimator::~RecursiveParameterEstimator()
{
}

RecursiveParameterEstimator::RecursiveParameterEstimator(const RecursiveParameterEstimator &another)
{
  this->make_copy(another);
}

void RecursiveParameterEstimator::operator=(const RecursiveParameterEstimator &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void RecursiveParameterEstimator::make_copy(const RecursiveParameterEstimator &another)
{
}
