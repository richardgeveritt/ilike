#include "gradient_estimator_output.h"

namespace ilike
{
GradientEstimatorOutput::GradientEstimatorOutput()
{
}

GradientEstimatorOutput::~GradientEstimatorOutput()
{
  
}

GradientEstimatorOutput::GradientEstimatorOutput(const GradientEstimatorOutput &another)
{
  this->make_copy(another);
}

void GradientEstimatorOutput::operator=(const GradientEstimatorOutput &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void GradientEstimatorOutput::make_copy(const GradientEstimatorOutput &another)
{
}
}
