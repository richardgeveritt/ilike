#include "exact_likelihood_estimator.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator(const NumericVector &data_in,
                                                   const EvaluateLogLikelihoodPtr &func_in)
:LikelihoodEstimator(data_in)
{
  this->func = func_in;
}

ExactLikelihoodEstimator::~ExactLikelihoodEstimator()
{

}

double ExactLikelihoodEstimator::estimate_log_likelihood(const NumericVector &inputs,
                                                         const List &auxiliary_variables) const
{
  return this->func(inputs,this->data);
}

List ExactLikelihoodEstimator::simulate_auxiliary_variables(const NumericVector &inputs) const
{
  return(0);
}

void ExactLikelihoodEstimator::setup_likelihood_estimator(const NumericMatrix &all_points,
                                                          const std::vector<List> &all_auxiliary_variables)
{

}
