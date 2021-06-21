#include "exact_likelihood_estimator.h"

ExactLikelihoodEstimator::ExactLikelihoodEstimator(const Data* observed_data_in,
                                                   const EvaluateLogLikelihoodPtr &func_in)
:LikelihoodEstimator(observed_data_in)
{
  this->func = func_in;
}

ExactLikelihoodEstimator::~ExactLikelihoodEstimator()
{

}

// double ExactLikelihoodEstimator::estimate_log_likelihood(const List &inputs,
//                                                          const List &auxiliary_variables) const
// {
//   return this->func(inputs,this->observed_data);
// }

LikelihoodEstimatorOutput* ExactLikelihoodEstimator::simulate(const Parameters &parameters) const
{
  return(new LikelihoodEstimatorOutput());
}

// void ExactLikelihoodEstimator::is_setup_likelihood_estimator(const std::vector<List> &all_points,
//                                                              const std::vector<List> &all_auxiliary_variables)
// {
//
// }
