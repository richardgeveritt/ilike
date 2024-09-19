#include "density_estimator_output.h"

namespace ilike
{
DensityEstimatorOutput::DensityEstimatorOutput()
{
  this->n = 0;
}

DensityEstimatorOutput::~DensityEstimatorOutput()
{
  
}

DensityEstimatorOutput::DensityEstimatorOutput(const DensityEstimatorOutput &another)
{
  this->make_copy(another);
}

void DensityEstimatorOutput::operator=(const DensityEstimatorOutput &another)
{
  if(this == &another)
    return;
  
  this->make_copy(another);
}

void DensityEstimatorOutput::make_copy(const DensityEstimatorOutput &another)
{
  this->n = another.n;
  //this->variables = another.variables;
}

void DensityEstimatorOutput::fit(const std::vector<Parameters> &points)
{
  arma::colvec normalised_log_weights(points.size(),arma::fill::value(-log(double(points.size()))));
  this->fit(points, normalised_log_weights);
}
}
