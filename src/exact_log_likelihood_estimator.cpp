#include "exact_log_likelihood_estimator.h"

ExactLogLikelihoodEstimator::ExactLogLikelihoodEstimator(const EvaluateLogLikelihoodPtr &func_in)
:LogLikelihoodEstimator()
{
  this->func = func_in;
}

ExactLogLikelihoodEstimator::~ExactLogLikelihoodEstimator()
{

}

double ExactLogLikelihoodEstimator::estimate_log_likelihood(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const
{
  return this->func(inputs,data);
}

List ExactLogLikelihoodEstimator::simulate_auxiliary_variables(const NumericVector &inputs,
                                                               const NumericVector &data) const
{
  return(0);
}

void ExactLogLikelihoodEstimator::setup_likelihood_estimator(const NumericMatrix &all_points,
                                                             const std::vector<List> &all_auxiliary_variables)
{

}
