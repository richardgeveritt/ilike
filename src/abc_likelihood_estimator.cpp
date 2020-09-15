#include "abc_likelihood_estimator.h"
#include "utils.h"

ABCLikelihoodEstimator::ABCLikelihoodEstimator(const NumericMatrix &data_in,
                                               const SimulateModelPtr &simulator_in,
                                               const double &number_of_likelihood_particles_in,
                                               const GetDataFromSimulationPtr &get_data_from_simulation_in,
                                               const EvaluateLogABCKernelPtr &evaluate_log_abc_kernel_in,
                                               const SummaryStatisticsPtr &summary_statistics_in,
                                               const double &abc_tolerance_in,
                                               const double &abc_desired_cess_in,
                                               const NumericVector &summary_statistics_scaling_in,
                                               const bool &adapt_abc_tolerance_to_cess_in,
                                               const bool &adapt_summary_statistics_scaling_in)
:LikelihoodEstimator(data_in)
{
  this->get_data_from_simulation = get_data_from_simulation_in;
  this->simulator = simulator_in;
  this->number_of_likelihood_particles = number_of_likelihood_particles_in;
  this->abc_likelihood = ABCLikelihood(evaluate_log_abc_kernel_in,
                                       summary_statistics_in,
                                       abc_tolerance_in,
                                       summary_statistics_scaling_in,
                                       data_in);
  this->adapt_summary_statistics_scaling = adapt_summary_statistics_scaling_in;
  this->abc_desired_cess = abc_desired_cess_in;
  this->adapt_abc_tolerance_to_cess = adapt_abc_tolerance_to_cess_in;
}

ABCLikelihoodEstimator::~ABCLikelihoodEstimator()
{

}

double ABCLikelihoodEstimator::estimate_log_likelihood(const NumericVector &inputs,
                                                       const List &auxiliary_variables) const
{
  std::vector<List> simulations = auxiliary_variables[0];
  std::vector<NumericMatrix> simulated_data;
  simulated_data.reserve(this->number_of_likelihood_particles);
  for (std::vector<List>::const_iterator i=simulations.begin(); i!=simulations.end(); ++i)
  {
    simulated_data.push_back(this->get_data_from_simulation(*i));
  }
  NumericVector abc_evaluations = this->abc_likelihood.evaluate_multiple(simulated_data);
  return log_sum_exp(abc_evaluations) - log(this->number_of_likelihood_particles);
}

List ABCLikelihoodEstimator::simulate_auxiliary_variables(const NumericVector &inputs) const
{
  std::vector<List> output;
  output.reserve(this->number_of_likelihood_particles);
  //NumericMatrix simulations(this->number_of_likelihood_particles,inputs.length());
  for (unsigned int i=0; i<this->number_of_likelihood_particles; ++i)
  {
    output.push_back(this->simulator(inputs, this->data));
  }

  return(List::create(output));
}

void ABCLikelihoodEstimator::is_setup_likelihood_estimator(const NumericMatrix &all_points,
                                                           const std::vector<List> &all_auxiliary_variables)
{

  if (this->adapt_summary_statistics_scaling)
  {
    unsigned int sum_stat_length = this->abc_likelihood.get_summary_data().length();

    // Use all_auxiliary_variables to find a good scaling for the summary statistics.
    NumericMatrix all_sumstats(all_auxiliary_variables.size()*this->number_of_likelihood_particles,
                               sum_stat_length);

    unsigned int counter = 0;

    std::vector< std::vector<NumericMatrix> > all_simulated_data;
    all_simulated_data.reserve(all_auxiliary_variables.size());
    for (std::vector<List>::const_iterator i=all_auxiliary_variables.begin(); i!=all_auxiliary_variables.end(); ++i)
    {
      std::vector<List> current_simulations = (*i)[0];
      std::vector<NumericMatrix> current_simulated_data;
      current_simulated_data.reserve(current_simulations.size());
      for (std::vector<List>::const_iterator j=current_simulations.begin(); j!=current_simulations.end(); ++j)
      {
        NumericMatrix a_simulation = this->get_data_from_simulation(*j);
        all_sumstats(counter,_) = this->abc_likelihood.summary_from_data(a_simulation);
        current_simulated_data.push_back(a_simulation);
        counter = counter + 1;
      }
      all_simulated_data.push_back(current_simulated_data);
    }

    NumericVector summary_statistics_scaling(sum_stat_length);
    for (unsigned int i=0; i<sum_stat_length; ++i)
    {
      summary_statistics_scaling[i] = 1.0/sd(all_sumstats(_,i));
    }

    this->abc_likelihood.set_summary_statistics_scaling(summary_statistics_scaling);
  }

  if (this->adapt_abc_tolerance_to_cess)
  {

    unsigned int n = all_auxiliary_variables.size();
    NumericVector current_log_weights(n);
    for (unsigned int i=0; i<n; ++i)
    {
      current_log_weights[i] = -log(double(n));
    }

    double abc_tolerance_local = epsilon_doubling(all_points,
                                                  all_auxiliary_variables,
                                                  current_log_weights);

    abc_tolerance_local = epsilon_bisection(abc_tolerance_local,
                                            all_points,
                                            all_auxiliary_variables,
                                            current_log_weights);

    this->abc_likelihood.set_abc_tolerance(abc_tolerance_local);
  }

}


double ABCLikelihoodEstimator::cess_score(const double &current_epsilon,
                                          const NumericMatrix &all_points,
                                          const std::vector<List> &all_auxiliary_variables,
                                          const NumericVector &current_log_weights)
{
  this->abc_likelihood.set_abc_tolerance(current_epsilon);
  unsigned int n = current_log_weights.length();
  NumericVector log_likelihoods(n);
  //unsigned int sum = 0;
  for (unsigned int i=0; i<n; ++i)
  {
    log_likelihoods[i] = this->estimate_log_likelihood(all_points(i,_),
                                                       all_auxiliary_variables[i]);

    //if (log_likelihoods[i] == R_NegInf)
    //  sum = sum + 1;
  }

  // if (sum==n)
  // {
  //   return(1.0 - this->abc_desired_cess);
  // }
  // else
  // {
  //
  // }
  return(cess(current_log_weights, log_likelihoods) - this->abc_desired_cess);
}


double ABCLikelihoodEstimator::epsilon_bisection(const double &current_epsilon,
                                                 const NumericMatrix &all_points,
                                                 const std::vector<List> &all_auxiliary_variables,
                                                 const NumericVector &current_log_weights)
{
  double current_bisect_epsilon = current_epsilon;
  double bisect_size = current_bisect_epsilon/2.0;
  double direction = -1;

  // double old_score = this->cess_score(current_bisect_epsilon,
  //                                     all_auxiliary_variables,
  //                                     current_log_weights);

  double new_bisect_epsilon;
  double the_score;

  for (unsigned int i=0; i<200; ++i)
  {
    new_bisect_epsilon = current_bisect_epsilon + direction*bisect_size;

    bisect_size = bisect_size/2.0;

    the_score = this->cess_score(new_bisect_epsilon,
                                 all_points,
                                 all_auxiliary_variables,
                                 current_log_weights);

    if (the_score >= 0) {
      direction = -1;
    } else {
      direction = 1;
    }

    current_bisect_epsilon = new_bisect_epsilon;

    if (the_score==0)
    {
      break;
    }

  }

  return(current_bisect_epsilon);
}


double ABCLikelihoodEstimator::epsilon_doubling(const NumericMatrix &all_points,
                                                const std::vector<List> &all_auxiliary_variables,
                                                const NumericVector &current_log_weights)
{
  double current_bisect_epsilon = 1;

  double old_score = this->cess_score(current_bisect_epsilon,
                                      all_points,
                                      all_auxiliary_variables,
                                      current_log_weights);

  if (old_score>=0)
    return(current_bisect_epsilon);

  double new_bisect_epsilon;
  double the_score;

  for (unsigned int i=0; i<100; ++i)
  {
    new_bisect_epsilon = 2.0*current_bisect_epsilon;

    the_score = cess_score(new_bisect_epsilon,
                           all_points,
                           all_auxiliary_variables,
                           current_log_weights);

    current_bisect_epsilon = new_bisect_epsilon;

    if (the_score>=0)
      return(current_bisect_epsilon);

  }

  throw Rcpp::exception("Doubling method cannot increase abc_tolerance sufficiently to reach desired CESS.");
}
