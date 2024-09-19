#include "direct_gradient_estimator.h"
#include "gradient_estimator_output.h"
#include "direct_gradient_estimator_output.h"

namespace ilike
{
DirectGradientEstimator::DirectGradientEstimator()
:GradientEstimator()
{
}

DirectGradientEstimator::~DirectGradientEstimator()
{
  
}

//Copy constructor for the DirectGradientEstimator class.
DirectGradientEstimator::DirectGradientEstimator(const DirectGradientEstimator &another)
:GradientEstimator(another)
{
  this->make_copy(another);
}

void DirectGradientEstimator::operator=(const DirectGradientEstimator &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  GradientEstimator::operator=(another);
  this->make_copy(another);
}

GradientEstimator* DirectGradientEstimator::duplicate() const
{
  return( new DirectGradientEstimator(*this));
}

void DirectGradientEstimator::make_copy(const DirectGradientEstimator &another)
{
}

GradientEstimatorOutput* DirectGradientEstimator::initialise()
{
  return new DirectGradientEstimatorOutput(this);
}

/*
 GradientEstimatorOutput* DirectGradientEstimator::generate_new_gradient_estimator_output(Particle &particle)
 {
 return new DirectGradientEstimatorOutput();
 }
 */

arma::mat DirectGradientEstimator::get_gradient_of_log(const std::string &variable,
                                                       const Index* index,
                                                       const Particle &particle) const
{
  return particle.direct_get_gradient_of_log(variable,
                                             index);
}

/*
 arma::mat DirectGradientEstimator::get_gradient_of_log(const std::string &variable,
 const Index* index,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 return particle.direct_get_gradient_of_log(variable,
 index,
 conditioned_on_parameters);
 }
 */

arma::mat DirectGradientEstimator::subsample_get_gradient_of_log(const std::string &variable,
                                                                 const Index* index,
                                                                 const Particle &particle) const
{
  return particle.direct_subsample_get_gradient_of_log(variable,
                                                       index);
}

/*
 arma::mat DirectGradientEstimator::subsample_get_gradient_of_log(const std::string &variable,
 const Index* index,
 Particle &particle,
 const Parameters &conditioned_on_parameters)
 {
 return particle.direct_subsample_get_gradient_of_log(variable,
 index,
 conditioned_on_parameters);
 }
 */
}
