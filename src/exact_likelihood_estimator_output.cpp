#include "exact_likelihood_estimator_output.h"
#include "exact_likelihood_estimator.h"
#include "utils.h"

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput()
  :LikelihoodEstimatorOutput()
{
}

ExactLikelihoodEstimatorOutput::ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in)
  :LikelihoodEstimatorOutput()
{
  this->estimator = estimator_in;
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
}

void ExactLikelihoodEstimatorOutput::continue_simulate(const Parameters &parameters)
{
}

void ExactLikelihoodEstimatorOutput::estimate(const Parameters &parameters)
{
  this->log_likelihood = this->estimator->func(parameters, *this->estimator->data);
}

void ExactLikelihoodEstimatorOutput::print(std::ostream &os) const
{

}
