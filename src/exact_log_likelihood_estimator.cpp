#include "exact_log_likelihood_estimator.h"

ExactLogLikelihoodEstimator::ExactLogLikelihoodEstimator(const EvaluateLogLikelihoodPtr &func_in)
:LogLikelihoodEstimator()
{
  this->func = func_in;
}

ExactLogLikelihoodEstimator::~ExactLogLikelihoodEstimator()
{

}

double ExactLogLikelihoodEstimator::operator()(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const
{
  return this->func(inputs,data);
}

List ExactLogLikelihoodEstimator::simulate_auxiliary_variables() const
{
  return(0);
}
