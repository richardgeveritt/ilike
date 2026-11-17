#include "density_estimator.h"

namespace ilike
{
DensityEstimator::DensityEstimator()
{
}

DensityEstimator::DensityEstimator(const std::vector<std::string> &variables_in)
{
  this->variables = variables_in;
}

DensityEstimator::~DensityEstimator()
{
  
}

DensityEstimator::DensityEstimator(const DensityEstimator &another)
{
  this->make_copy(another);
}

void DensityEstimator::operator=(const DensityEstimator &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void DensityEstimator::make_copy(const DensityEstimator &another)
{
  this->variables = another.variables;
}

std::vector<std::string> DensityEstimator::get_variables() const
{
  return this->variables;
}

}
