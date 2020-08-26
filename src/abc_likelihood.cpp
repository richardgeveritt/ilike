#include "abc_likelihood.h"

ABCLikelihood::ABCLikelihood()
{
}

ABCLikelihood::ABCLikelihood(const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                             const SummaryStatisticPtr &summary_statistic_in,
                             const NumericVector &abc_tolerances_in,
                             const NumericVector &summary_statistic_scaling_in,
                             const NumericVector &data_in)
{
  this->evaluate_log_abc_kernel = evaluate_log_abc_kernel_in;
  this->summary_statistic = summary_statistic_in;
  this->abc_tolerances = abc_tolerances_in;
  this->summary_data = this->summary_statistic(data_in);
}

ABCLikelihood::ABCLikelihood(const ABCLikelihood &another)
{
  this->MakeCopy(another);
}

ABCLikelihood::~ABCLikelihood()
{

}

void ABCLikelihood::operator=(const ABCLikelihood &another)
{
  if(this == &another){ //if a==a
    return;
  }

  this->MakeCopy(another);
}

void ABCLikelihood::MakeCopy(const ABCLikelihood &another)
{
  this->abc_tolerances = another.abc_tolerances;
  this->evaluate_log_abc_kernel = another.evaluate_log_abc_kernel;
  this->summary_statistic = another.summary_statistic;
  this->summary_statistic_scaling = another.summary_statistic_scaling;
  this->summary_data = another.summary_data;
}

double ABCLikelihood::evaluate(const NumericVector &simulated_data) const
{
  return this->evaluate_log_abc_kernel(this->summary_statistic(simulated_data), this->summary_data, this->abc_tolerances);
}

NumericVector ABCLikelihood::get_summary_data() const
{
  return this->summary_data;
}

SummaryStatisticPtr ABCLikelihood::get_summary_statistic() const
{
  return this->summary_statistic;
}

void ABCLikelihood::set_abc_tolerances(const NumericVector &abc_tolerances_in)
{
  this->abc_tolerances = abc_tolerances_in;
}

void ABCLikelihood::set_summary_statistic_scaling(const NumericVector &summary_statistic_scaling_in)
{
  this->summary_statistic_scaling = summary_statistic_scaling_in;
}
