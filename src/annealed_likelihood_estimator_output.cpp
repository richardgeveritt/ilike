#include "annealed_likelihood_estimator_output.h"
#include "annealed_likelihood_estimator.h"
#include "utils.h"

AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
}

AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput(AnnealedLikelihoodEstimator* estimator_in,
     LikelihoodEstimatorOutput* estimator_output_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->estimator_output = estimator_output_in;
}

AnnealedLikelihoodEstimatorOutput::~AnnealedLikelihoodEstimatorOutput()
{
  if (this->estimator_output!=NULL)
    delete this->estimator_output;
}

//Copy constructor for the AnnealedLikelihoodEstimatorOutput class.
AnnealedLikelihoodEstimatorOutput::AnnealedLikelihoodEstimatorOutput(const AnnealedLikelihoodEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void AnnealedLikelihoodEstimatorOutput::operator=(const AnnealedLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->estimator_output!=NULL)
    delete this->estimator_output;

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* AnnealedLikelihoodEstimatorOutput::duplicate(void)const
{
  return( new AnnealedLikelihoodEstimatorOutput(*this));
}

void AnnealedLikelihoodEstimatorOutput::make_copy(const AnnealedLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  if (another.estimator_output!=NULL)
    this->estimator_output = another.estimator_output->duplicate();
  else
    this->estimator_output = NULL;
}

void AnnealedLikelihoodEstimatorOutput::simulate()
{
  this->estimator_output->simulate();
}

void AnnealedLikelihoodEstimatorOutput::simulate(const Parameters &parameters)
{
  this->estimator_output->simulate(parameters);
}

void AnnealedLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator_output->evaluate(parameters);
  }
  else
  {
    this->estimator_output->evaluate_smcfixed_part(parameters);
  }
}

void AnnealedLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator_output->evaluate_smcadaptive_part_given_smcfixed(parameters);
  }
  
  if (this->estimator->use_constant)
    this->log_likelihood = this->estimator->constant_power * this->estimator_output->log_likelihood;
  else
  {
    this->log_likelihood = this->estimator->function_power(parameters) * this->estimator_output->log_likelihood;
  }
}

void AnnealedLikelihoodEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
  this->estimator_output->subsample_simulate(parameters);
}

void AnnealedLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
  {
    this->estimator_output->subsample_evaluate(parameters);
  }
  else
  {
    this->estimator_output->subsample_evaluate_smcfixed_part(parameters);
  }
}

void AnnealedLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (!this->estimator->smcfixed_flag)
  {
    this->estimator_output->subsample_evaluate_smcadaptive_part_given_smcfixed(parameters);
  }
  
  if (this->estimator->use_constant)
    this->log_likelihood = this->estimator->constant_power * this->estimator_output->subsample_log_likelihood;
  else
  {
    this->log_likelihood = this->estimator->function_power(parameters) * this->estimator_output->subsample_log_likelihood;
  }
}

LikelihoodEstimator* AnnealedLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat AnnealedLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                                 const Parameters &x)
{
  arma::mat gradient;
  if (this->estimator->use_constant)
    gradient = this->estimator->constant_power * this->estimator_output->get_gradient_of_log(variable,
                                                                                             x);
  else
  {
    gradient = this->estimator->function_power(x) * this->estimator_output->get_gradient_of_log(variable,
                                                                                                x);
  }
  return gradient;
}

arma::mat AnnealedLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                 const Parameters &x)
{
  arma::mat gradient;
  if (this->estimator->use_constant)
    gradient = this->estimator->constant_power * this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                                             x);
  else
  {
    gradient = this->estimator->function_power(x) * this->estimator_output->subsample_get_gradient_of_log(variable,
                                                                                                x);
  }
  return gradient;
}

void AnnealedLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}
