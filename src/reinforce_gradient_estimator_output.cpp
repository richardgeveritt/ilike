#include "reinforce_gradient_estimator_output.h"
#include "reinforce_gradient_estimator.h"

ReinforceInfo::ReinforceInfo()
{
}

ReinforceInfo::~ReinforceInfo()
{
}

ReinforceInfo::ReinforceInfo(const ReinforceInfo &another)
{
  this->make_copy(another);
}

void ReinforceInfo::operator=(const ReinforceInfo &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->make_copy(another);
}

void ReinforceInfo::make_copy(const ReinforceInfo &another)
{
  this->gradient = another.gradient;
}

ReinforceGradientEstimatorOutput::ReinforceGradientEstimatorOutput()
  :GradientEstimatorOutput()
{
  this->estimator = NULL;
}

ReinforceGradientEstimatorOutput::ReinforceGradientEstimatorOutput(ReinforceGradientEstimator* estimator_in)
{
  this->estimator = estimator_in;
}

ReinforceGradientEstimatorOutput::~ReinforceGradientEstimatorOutput()
{
  
}

// Copy constructor for the ReinforceGradientEstimatorOutput class.
ReinforceGradientEstimatorOutput::ReinforceGradientEstimatorOutput(const ReinforceGradientEstimatorOutput &another)
  :GradientEstimatorOutput(another)
{
  this->make_copy(another);
}

void ReinforceGradientEstimatorOutput::operator=(const ReinforceGradientEstimatorOutput &another)
{
  if(this == &another)
  {
    return;
  }
  
  GradientEstimatorOutput::operator=(another);
  this->make_copy(another);
}

GradientEstimatorOutput* ReinforceGradientEstimatorOutput::duplicate() const
{
  return( new ReinforceGradientEstimatorOutput(*this));
}

void ReinforceGradientEstimatorOutput::make_copy(const ReinforceGradientEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->infos = another.infos;
  this->subsample_infos = another.infos;
}

arma::mat ReinforceGradientEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                                const Index* index,
                                                                const Particle &particle)
{
  auto found = this->infos.find(variable);
  
  if (found != this->infos.end())
  {
    return found->second.gradient;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    ReinforceInfo info = ReinforceInfo();
    info.gradient = this->estimator->get_gradient_of_log(variable,
                                                         this->auxiliary_variables.find(variable)->second,
                                                         index,
                                                         particle);
    this->infos[variable] = info;
    return info.gradient;
  }
}

/*
arma::mat ReinforceGradientEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                                const Index* index,
                                                                Particle &particle,
                                                                const Parameters &conditioned_on_parameters)
{
  auto found = this->infos.find(variable);
  
  if (found != this->infos.end())
  {
    return found->second.gradient;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    ReinforceInfo info = ReinforceInfo();
    info.gradient = this->estimator->get_gradient_of_log(variable,index,particle,conditioned_on_parameters);
    this->infos[variable] = info;
    return info.gradient;
  }
}
*/

arma::mat ReinforceGradientEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                          const Index* index,
                                                                          const Particle &particle)
{
  auto found = this->infos.find(variable);
  
  if (found != this->infos.end())
  {
    return found->second.gradient;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    ReinforceInfo info = ReinforceInfo();
    info.gradient = this->estimator->subsample_get_gradient_of_log(variable,
                                                                   this->auxiliary_variables.find(variable)->second,
                                                                   index,
                                                                   particle);
    this->infos[variable] = info;
    return info.gradient;
  }
}

void ReinforceGradientEstimatorOutput::simulate_auxiliary_variables()
{
  this->auxiliary_variables = this->estimator->simulate_auxiliary_variables();
}

/*
arma::mat ReinforceGradientEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                          const Index* index,
                                                                          Particle &particle,
                                                                          const Parameters &conditioned_on_parameters)
{
  auto found = this->infos.find(variable);
  
  if (found != this->infos.end())
  {
    return found->second.gradient;
  }
  else
  {
    // if proposal is not found, regenerate all of the variables needed to estimate the gradient
    ReinforceInfo info = ReinforceInfo();
    info.gradient = this->estimator->subsample_get_gradient_of_log(variable,index,particle,conditioned_on_parameters);
    this->infos[variable] = info;
    return info.gradient;
  }
}
*/
