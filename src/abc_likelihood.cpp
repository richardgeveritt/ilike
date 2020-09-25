#include "abc_likelihood.h"
#include "utils.h"


ABCLikelihood::ABCLikelihood()
{
}

ABCLikelihood::ABCLikelihood(const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                             const SummaryStatisticsPtr &summary_statistics_in,
                             const double &abc_tolerance_in,
                             const arma::colvec &summary_statistics_scaling_in,
                             const List &observed_data_in)
{
  this->evaluate_log_abc_kernel = evaluate_log_abc_kernel_in;
  this->summary_statistics = summary_statistics_in;
  this->summary_statistics_scaling = summary_statistics_scaling_in;
  this->abc_tolerance = abc_tolerance_in;
  this->summary_data = this->summary_statistics(observed_data_in);
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
  this->abc_tolerance = another.abc_tolerance;
  this->evaluate_log_abc_kernel = another.evaluate_log_abc_kernel;
  this->summary_statistics = another.summary_statistics;
  this->summary_statistics_scaling = another.summary_statistics_scaling;
  this->summary_data = another.summary_data;
}

double ABCLikelihood::evaluate(const List &simulated_data) const
{
  arma::colvec simulated_summary_stats = this->summary_statistics_scaling % this->summary_statistics(simulated_data);
  arma::colvec observed_summary_stats = this->summary_statistics_scaling % this->summary_data;
  double p = 2.0;
  unsigned int n = simulated_summary_stats.size();
  double Lp_log_distance = (1.0/p)*log_sum_exp(p*log(abs(simulated_summary_stats-observed_summary_stats)));
  double Lp_ball_log_volume = double(n)*(log(2.0) + lgamma((1.0/p)+1.0) + log(abc_tolerance)) - lgamma((double(n)/p)+1.0);
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume);
  else
    return(R_NegInf);

  // return this->evaluate_log_abc_kernel(this->summary_statistics_scaling % this->summary_statistics(simulated),
  //                                      this->summary_statistics_scaling % this->summary_data,
  //                                      this->abc_tolerance);
}

arma::colvec ABCLikelihood::evaluate_multiple(const std::vector<List> &simulations) const
{
  unsigned int n = simulations.size();
  arma::colvec log_likelihoods(n);
  unsigned int counter = 0;

  for (std::vector<List>::const_iterator i=simulations.begin(); i!=simulations.end(); ++i)
  {
    log_likelihoods[counter] = this->evaluate(*i);
    counter = counter + 1;
  }

  return log_likelihoods;
}

arma::colvec ABCLikelihood::summary_from_data(const List &simulation) const
{
  return(this->summary_statistics(simulation));
}

arma::colvec ABCLikelihood::get_summary_data() const
{
  return this->summary_data;
}

SummaryStatisticsPtr ABCLikelihood::get_summary_statistics() const
{
  return this->summary_statistics;
}

double ABCLikelihood::get_abc_tolerance() const
{
  return this->abc_tolerance;
}

void ABCLikelihood::set_abc_tolerance(const double &abc_tolerance_in)
{
  this->abc_tolerance = abc_tolerance_in;
}

void ABCLikelihood::set_summary_statistics_scaling(const arma::colvec &summary_statistics_scaling_in)
{
  this->summary_statistics_scaling = summary_statistics_scaling_in;
}
