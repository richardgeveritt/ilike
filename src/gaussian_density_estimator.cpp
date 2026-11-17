#include "gaussian_density_estimator.h"
#include "sample_mean_vector_parameter_estimator.h"
#include "sample_covariance_matrix_parameter_estimator.h"
#include "distributions.h"
#include "gaussian_density_estimator_output.h"

namespace ilike
{
GaussianDensityEstimator::GaussianDensityEstimator()
:DensityEstimator()
{
  this->unbiased = true;
}

GaussianDensityEstimator::GaussianDensityEstimator(const std::vector<std::string> &variables_in)
:DensityEstimator(variables_in)
{
}

GaussianDensityEstimator::GaussianDensityEstimator(const std::vector<std::string> &variables_in,
                                                   bool unbiased_in)
:DensityEstimator(variables_in)
{
  this->unbiased = unbiased_in;
  
}

GaussianDensityEstimator::~GaussianDensityEstimator()
{
  
}

//Copy constructor for the GaussianDensityEstimator class.
GaussianDensityEstimator::GaussianDensityEstimator(const GaussianDensityEstimator &another)
:DensityEstimator(another)
{
  this->make_copy(another);
}

void GaussianDensityEstimator::operator=(const GaussianDensityEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }

  DensityEstimator::operator=(another);
  this->make_copy(another);
}

DensityEstimator* GaussianDensityEstimator::duplicate() const
{
  return( new GaussianDensityEstimator(*this));
}

void GaussianDensityEstimator::make_copy(const GaussianDensityEstimator &another)
{
  this->unbiased = another.unbiased;
  
}

DensityEstimatorOutput* GaussianDensityEstimator::initialise()
{
  return new GaussianDensityEstimatorOutput(this);
}

bool GaussianDensityEstimator::get_unbiased() const
{
  return this->unbiased;
}

}
