#include "direct_gradient_estimator_output.h"
#include "direct_gradient_estimator.h"
#include "proposal_kernel.h"
#include "transform.h"

DirectGradientEstimatorOutput::DirectGradientEstimatorOutput()
  :GradientEstimatorOutput()
{
  this->estimator = NULL;
}

DirectGradientEstimatorOutput::DirectGradientEstimatorOutput(DirectGradientEstimator* estimator_in)
{
  this->estimator = estimator_in;
}

DirectGradientEstimatorOutput::~DirectGradientEstimatorOutput()
{
  
}

//Copy constructor for the DirectGradientEstimatorOutput class.
DirectGradientEstimatorOutput::DirectGradientEstimatorOutput(const DirectGradientEstimatorOutput &another)
  :GradientEstimatorOutput(another)
{
  this->make_copy(another);
}

void DirectGradientEstimatorOutput::operator=(const DirectGradientEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  GradientEstimatorOutput::operator=(another);
  this->make_copy(another);
}

GradientEstimatorOutput* DirectGradientEstimatorOutput::duplicate() const
{
  return( new DirectGradientEstimatorOutput(*this));
}

void DirectGradientEstimatorOutput::make_copy(const DirectGradientEstimatorOutput &another)
{
  this->gradients = another.gradients;
  this->subsample_gradients = another.subsample_gradients;
  this->estimator = another.estimator;
}

arma::mat DirectGradientEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                             const Index* index,
                                                             const Particle &particle)
{
  auto found = this->gradients.find(variable);
  
  if (found != this->gradients.end())
  {
    return found->second;
  }
  else
  {

    arma::mat current_gradient = particle.direct_get_gradient_of_log(variable,
                                                                     index);
    if (this->estimator->proposal->transform!=NULL)
    {
      current_gradient = arma::reshape(arma::vectorise(current_gradient)*this->estimator->proposal->transform->inverse_jacobian(particle.get_transformed_parameters(this->estimator->proposal)),current_gradient.n_rows,current_gradient.n_cols);
    }
    this->gradients[variable] = current_gradient;
    return current_gradient;
  }
}

/*
arma::mat DirectGradientEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                             const Index* index,
                                                             Particle &particle,
                                                             const Parameters &conditioned_on_parameters_in)
{
  auto found = this->gradients.find(variable);
  
  if (found != this->gradients.end())
  {
    return found->second;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    arma::mat current_gradient = particle.direct_get_gradient_of_log(variable,
                                                                     index,
                                                                     conditioned_on_parameters_in);
    if (this->estimator->proposal->transform!=NULL)
    {
      current_gradient = arma::reshape(arma::vectorise(current_gradient)*this->estimator->proposal->transform->inverse_jacobian(*particle.move_parameters),current_gradient.n_rows,current_gradient.n_cols);
    }
    this->gradients[variable] = current_gradient;
    return current_gradient;
  }
}
*/

arma::mat DirectGradientEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                       const Index* index,
                                                                       const Particle &particle)
{
  auto found = this->gradients.find(variable);
  
  if (found != this->gradients.end())
  {
    return found->second;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    arma::mat current_gradient = particle.direct_subsample_get_gradient_of_log(variable,
                                                                               index);
    if (this->estimator->proposal->transform!=NULL)
    {
      current_gradient = arma::reshape(arma::vectorise(current_gradient)*this->estimator->proposal->transform->inverse_jacobian(particle.get_transformed_parameters(this->estimator->proposal)),current_gradient.n_rows,current_gradient.n_cols);
    }
    this->gradients[variable] = current_gradient;
    return current_gradient;
  }
}

void DirectGradientEstimatorOutput::simulate_auxiliary_variables()
{
  
}

/*
arma::mat DirectGradientEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                       const Index* index,
                                                                       Particle &particle,
                                                                       const Parameters &conditioned_on_parameters_in)
{
  auto found = this->gradients.find(variable);
  
  if (found != this->gradients.end())
  {
    return found->second;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    arma::mat current_gradient = particle.direct_subsample_get_gradient_of_log(variable,
                                                                               index,
                                                                               conditioned_on_parameters_in);
    if (this->estimator->proposal->transform!=NULL)
    {
      current_gradient = arma::reshape(arma::vectorise(current_gradient)*this->estimator->proposal->transform->inverse_jacobian(*particle.move_parameters),current_gradient.n_rows,current_gradient.n_cols);
    }
    this->gradients[variable] = current_gradient;
    return current_gradient;
  }
}
*/
