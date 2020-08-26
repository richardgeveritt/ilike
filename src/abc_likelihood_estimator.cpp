#include "abc_likelihood_estimator.h"
#include "utils.h"

ABCLikelihoodEstimator::ABCLikelihoodEstimator(const NumericVector &data_in,
                                               const SimulateModelPtr &simulator_in,
                                               const double &number_of_simulations_in,
                                               const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                                               const SummaryStatisticPtr &summary_statistic_in,
                                               const NumericVector &tolerances_in,
                                               const NumericVector &summary_statistic_scaling_in)
:LikelihoodEstimator(data_in)
{
  this->simulator = simulator_in;
  this->number_of_simulations = number_of_simulations_in;
  this->abc_likelihood = ABCLikelihood(evaluate_log_abc_kernel_in,
                                       summary_statistic_in,
                                       tolerances_in,
                                       summary_statistic_scaling_in,
                                       data_in);
}

ABCLikelihoodEstimator::~ABCLikelihoodEstimator()
{

}

double ABCLikelihoodEstimator::estimate_log_likelihood(const NumericVector &inputs,
                                                       const List &auxiliary_variables) const
{
  NumericMatrix aux_variables = auxiliary_variables[0];
  NumericVector abc_evaluations(this->number_of_simulations);
  for (unsigned int i=0; i<this->number_of_simulations; ++i)
  {
    abc_evaluations[i] = this->abc_likelihood.evaluate(aux_variables(i,_));
  }
  return log_sum_exp(abc_evaluations) - log(this->number_of_simulations);
}

List ABCLikelihoodEstimator::simulate_auxiliary_variables(const NumericVector &inputs) const
{
  NumericMatrix simulations(this->number_of_simulations,inputs.length());
  for (unsigned int i=0; i<this->number_of_simulations; ++i)
  {
    simulations(i,_) = this->simulator(inputs,this->data);
  }

  return(List::create(simulations));
}

void ABCLikelihoodEstimator::setup_likelihood_estimator(const NumericMatrix &all_points,
                                                        const std::vector<List> &all_auxiliary_variables)
{
  if (all_auxiliary_variables.size()==0)
    return;

  unsigned int sum_stat_length = this->abc_likelihood.get_summary_data().length();

  // Use all_auxiliary_variables to find a good scaling for the summary statistics.
  NumericMatrix all_sumstats(all_auxiliary_variables.size()*this->number_of_simulations,
                             sum_stat_length);

  unsigned int counter = 0;

  for (unsigned int i=0; i<all_auxiliary_variables.size(); ++i)
  {
    NumericMatrix current_aux_variables = all_auxiliary_variables[i][0];
    for (unsigned int j=0; this->number_of_simulations; ++j)
    {
      all_sumstats(counter,_) = this->abc_likelihood.get_summary_statistic()(current_aux_variables(j,_));
      counter = counter + 1;
    }
  }

  NumericVector summary_statistic_scaling(sum_stat_length);
  for (unsigned int i=0; i<sum_stat_length; ++i)
  {
    summary_statistic_scaling[i] = sd(all_sumstats(_,i));
  }

  this->abc_likelihood.set_summary_statistic_scaling(summary_statistic_scaling);
}
