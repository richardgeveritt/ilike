#include "exact_likelihood_estimator.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator(const List &observed_data_in,
                                                   const EvaluateLogLikelihoodPtr &func_in)
:LikelihoodEstimator(observed_data_in)
{
  this->func = func_in;
}

ExactLikelihoodEstimator::~ExactLikelihoodEstimator()
{

}

double ExactLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
                                                         const List &auxiliary_variables) const
{
  return this->func(inputs,this->observed_data);
}

List ExactLikelihoodEstimator::simulate_auxiliary_variables(const List &inputs) const
{
  return(0);
}

void ExactLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
                                                             const std::vector<List> &all_auxiliary_variables)
{

}
