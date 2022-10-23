#include "exact_likelihood_estimator_output.h"
#include "exact_likelihood_estimator.h"
#include "utils.h"

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
  this->log_likelihood_smcfixed_part = 0.0;
  this->subsample_log_likelihood_smcfixed_part = 0.0;
}

ExactLikelihoodEstimatorOutput::~ExactLikelihoodEstimatorOutput()
{

}

//Copy constructor for the ExactLikelihoodEstimatorOutput class.
ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput(const ExactLikelihoodEstimatorOutput &another)
  :LikelihoodEstimatorOutput(another)
{
  this->make_copy(another);
}

void ExactLikelihoodEstimatorOutput::operator=(const ExactLikelihoodEstimatorOutput &another)
{
  if(this == &another){ //if a==a
    return;
  }

  LikelihoodEstimatorOutput::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimatorOutput* ExactLikelihoodEstimatorOutput::duplicate(void)const
{
  return( new ExactLikelihoodEstimatorOutput(*this));
}

void ExactLikelihoodEstimatorOutput::make_copy(const ExactLikelihoodEstimatorOutput &another)
{
  this->estimator = another.estimator;
  this->log_likelihood_smcfixed_part = another.log_likelihood_smcfixed_part;
  this->subsample_log_likelihood_smcfixed_part = another.subsample_log_likelihood_smcfixed_part;
}

void ExactLikelihoodEstimatorOutput::simulate()
{
}

void ExactLikelihoodEstimatorOutput::simulate(const Parameters &parameters)
{
}

void ExactLikelihoodEstimatorOutput::evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood_smcfixed_part = this->estimator->evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->log_likelihood = this->log_likelihood_smcfixed_part;
  else
    this->log_likelihood = this->estimator->evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::subsample_simulate(const Parameters &parameters)
{
}

void ExactLikelihoodEstimatorOutput::subsample_evaluate_smcfixed_part(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood_smcfixed_part = this->estimator->subsample_evaluate(parameters);
}

void ExactLikelihoodEstimatorOutput::subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)
{
  if (this->estimator->smcfixed_flag)
    this->subsample_log_likelihood = this->subsample_log_likelihood_smcfixed_part;
  else
    this->subsample_log_likelihood = this->estimator->subsample_evaluate(parameters);
}

LikelihoodEstimator* ExactLikelihoodEstimatorOutput::get_likelihood_estimator() const
{
  return this->estimator;
}

arma::mat ExactLikelihoodEstimatorOutput::get_gradient_of_log(const std::string &variable,
                                                              const Parameters &x)
{
  return this->estimator->evaluate_gradient(variable,x);
}

arma::mat ExactLikelihoodEstimatorOutput::subsample_get_gradient_of_log(const std::string &variable,
                                                                        const Parameters &x)
{
  return this->estimator->subsample_evaluate_gradient(variable,x);
}

void ExactLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}
