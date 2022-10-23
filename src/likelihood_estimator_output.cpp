#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "data.h"

LikelihoodEstimatorOutput::LikelihoodEstimatorOutput()
{
  this->log_likelihood = 0.0;
  this->subsample_log_likelihood = 0.0;
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
  this->subsample_log_likelihood = another.subsample_log_likelihood;
}

double LikelihoodEstimatorOutput::evaluate(const Parameters &parameters)
{
  this->evaluate_smcfixed_part(parameters);
  this->evaluate_smcadaptive_part_given_smcfixed(parameters);
  return this->log_likelihood;
}

double LikelihoodEstimatorOutput::subsample_evaluate(const Parameters &parameters)
{
  this->subsample_evaluate_smcfixed_part(parameters);
  this->subsample_evaluate_smcadaptive_part_given_smcfixed(parameters);
  return this->subsample_log_likelihood;
}

void LikelihoodEstimatorOutput::change_data()
{
  this->get_likelihood_estimator()->change_data();
}

void LikelihoodEstimatorOutput::change_data(Data* new_data)
{
  this->get_likelihood_estimator()->change_data(new_data);
}

std::ostream& operator<<(std::ostream& os, const LikelihoodEstimatorOutput &output)
{
  output.print(os);
  return os;
}

void LikelihoodEstimatorOutput::print(std::ostream &os) const
{
}
