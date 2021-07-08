#include "likelihood_estimator_output.h"

LikelihoodEstimatorOutput::LikelihoodEstimatorOutput()
{
}

LikelihoodEstimatorOutput::~LikelihoodEstimatorOutput()
{

}

// void LikelihoodEstimatorOutput::simulate(const Parameters &parameters)
// {
//
// }

LikelihoodEstimatorOutput::LikelihoodEstimatorOutput(const LikelihoodEstimatorOutput &another)
{
  this->make_copy(another);
}

void LikelihoodEstimatorOutput::operator=(const LikelihoodEstimatorOutput &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void LikelihoodEstimatorOutput::make_copy(const LikelihoodEstimatorOutput &another)
{
  this->log_likelihood = another.log_likelihood;
}
