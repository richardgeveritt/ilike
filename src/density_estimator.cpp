#include "density_estimator.h"

DensityEstimator::DensityEstimator()
{
  //this->n = 0;
}

DensityEstimator::DensityEstimator(const std::vector<std::string> &variables_in)
{
  //this->n = 0;
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
  //this->n = another.n;
  this->variables = another.variables;
}

std::vector<std::string> DensityEstimator::get_variables() const
{
  return this->variables;
}

/*
void DensityEstimator::fit(const std::vector<Parameters> &points)
{
  arma::colvec normalised_log_weights(points.size(),arma::fill::value(-log(double(points.size()))));
  this->fit(points, normalised_log_weights);
}
*/
