Lp_uniform_evaluate_log_abc_kernel = function(simulated_summary_stats,
                                              observed_summary_stats,
                                              abc_tolerance,
                                              p)
{
  n = length(simulated_summary_stats)
  Lp_log_distance = (1/p)*log_sum_exp(p*log(abs(simulated_summary_stats-observed_summary_stats)))
  Lp_ball_log_volume = n*(log(2) + lgamma((1/p)+1) + log(abc_tolerance)) - lgamma((n/p)+1)
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume)
  else
    return(-Inf)
}


Linf_uniform_evaluate_log_abc_kernel = function(simulated_summary_stats,
                                                observed_summary_stats,
                                                abc_tolerance)
{
  n = length(simulated_summary_stats)
  Lp_log_distance = log(max(abs(simulated_summary_stats-observed_summary_stats)))
  Lp_ball_log_volume = n*(log(2) + log(abc_tolerance))
  if (Lp_log_distance<=log(abc_tolerance))
    return(-Lp_ball_log_volume)
  else
    return(-Inf)
}


gaussian_evaluate_log_abc_kernel = function(simulated_summary_stats,
                                            observed_summary_stats,
                                            abc_tolerance)
{
  log_probs = apply(simulated_summary_stats-observed_summary_stats, 1, dnorm, 0, abc_tolerance, TRUE)
  return(sum(log_probs))
}


evaluate_abc_likelihood = function(simulated_data,
                                   evaluate_log_abc_kernel,
                                   summary_statistics,
                                   summary_data,
                                   abc_tolerance,
                                   summary_statistics_scaling)
{
  return(evaluate_log_abc_kernel(summary_statistics_scaling*summary_statistics(simulated_data),
                                 summary_statistics_scaling*summary_data,
                                 abc_tolerance))
}

multiple_evaluate_abc_likelihood = function(simulations,
                                            evaluate_log_abc_kernel,
                                            summary_statistics,
                                            summary_data,
                                            abc_tolerance,
                                            summary_statistics_scaling)
{
  if (likelihood_use_future==TRUE)
  {
    return(unlist(future.apply::future_lapply(simulations,
                                              FUN = function(i){
                                                evaluate_abc_likelihood(i,
                                                                        evaluate_log_abc_kernel,
                                                                        summary_statistics,
                                                                        summary_data,
                                                                        abc_tolerance,
                                                                        summary_statistics_scaling)})))
  }
  else
  {
    return(unlist(lapply(simulations, FUN = function(i){
      evaluate_abc_likelihood(i,
                              evaluate_log_abc_kernel,
                              summary_statistics,
                              summary_data,
                              abc_tolerance,
                              summary_statistics_scaling)})))
  }
}

abc_simulate_auxiliary_variables = function(inputs,
                                            number_of_likelihood_particles,
                                            simulator,
                                            observed_data)
{
  repeated_inputs = lapply(1:number_of_likelihood_particles, FUN = function(i){inputs})

  if (likelihood_use_future==TRUE)
  {
    simulations = future.apply::future_lapply(repeated_inputs,
                              FUN = function(i){simulator(i, observed_data)},
                              future.seed = TRUE)
  }
  else
  {
    simulations = lapply(repeated_inputs,
                         FUN = function(i){simulator(i, observed_data)})
  }

  return(list(simulations))
}

abc_estimate_log_likelihood = function(auxiliary_variables,
                                       evaluate_log_abc_kernel,
                                       summary_statistics,
                                       summary_data,
                                       abc_tolerance,
                                       summary_statistics_scaling)
{
  simulations = auxiliary_variables[[1]]
  #simulated_data = lapply(simulations, function(s){return(s["data"])})
  abc_evaluations = multiple_evaluate_abc_likelihood(simulations,
                                                     evaluate_log_abc_kernel,
                                                     summary_statistics,
                                                     summary_data,
                                                     abc_tolerance,
                                                     summary_statistics_scaling)

  return(log_sum_exp(abc_evaluations) - log(length(abc_evaluations)))
}

abc_is_setup_likelihood_estimator = function(all_points,
                                             all_auxiliary_variables,
                                             number_of_likelihood_particles,
                                             evaluate_log_abc_kernel,
                                             summary_statistics,
                                             summary_data,
                                             abc_tolerance,
                                             abc_desired_cess,
                                             summary_statistics_scaling,
                                             adapt_abc_tolerance_to_cess,
                                             adapt_summary_statistics_scaling)
{
  # Resulting function is the standard ABC estimator: it will take some simulated data variables and find the average of an ABC kernel at these values.

  if (adapt_summary_statistics_scaling)
  {
    sum_stat_length = length(summary_data)

    if (length(all_auxiliary_variables)==0)
    {
      summary_statistics_scaling = matrix(1,sum_stat_length)
    }
    else
    {

      #data_length = ncol(all_auxiliary_variables[[1]])

      # Use all_auxiliary_variables to find a good scaling for the summary statistics.

      list_of_list_of_sumstats = lapply(all_auxiliary_variables, FUN = function(aux_vars){
         return(lapply(aux_vars[[1]], FUN = function(a_simulation){
          summary_statistics(a_simulation)
        }))
      })

      all_summary_statistics = t(matrix(unlist(list_of_list_of_sumstats),sum_stat_length,number_of_likelihood_particles*length(all_auxiliary_variables)))

      #sumstats_vec = sapply(sumstats,c)

      summary_statistics_scaling = 1/apply(all_summary_statistics, 2, sd, na.rm=TRUE)

    }
  }

  if (adapt_abc_tolerance_to_cess==TRUE)
  {

    all_abc_estimate_log_likelihoods_tol = function(tol)
    {
      return(unlist(lapply(1:length(all_auxiliary_variables),function(i){
        abc_estimate_log_likelihood(all_auxiliary_variables[[i]], evaluate_log_abc_kernel, summary_statistics, summary_data, tol, summary_statistics_scaling) })))
    }

    current_log_weights = matrix(-log(length(all_auxiliary_variables)), length(all_auxiliary_variables))

    # Use doubling to find abc_tolerance that gives ESS above desired ESS.
    abc_tolerance = epsilon_doubling(all_auxiliary_variables,
                                     abc_desired_cess,
                                     current_log_weights,
                                     all_abc_estimate_log_likelihoods_tol)

    # Use bisection to find abc_tolerance that gives desired ESS.
    abc_tolerance = epsilon_bisection(abc_tolerance,
                                      all_auxiliary_variables,
                                      abc_desired_cess,
                                      current_log_weights,
                                      all_abc_estimate_log_likelihoods_tol)

  }

  return(function(point,auxiliary_variables){return(abc_estimate_log_likelihood(auxiliary_variables,
                                                                                evaluate_log_abc_kernel,
                                                                                summary_statistics,
                                                                                summary_data,
                                                                                abc_tolerance,
                                                                                summary_statistics_scaling))})

}


cess_score <- function(current_epsilon,
                       all_auxiliary_variables,
                       abc_desired_cess,
                       current_log_weights,
                       all_abc_estimate_log_likelihoods_tol)
{
  log_likelihoods = all_abc_estimate_log_likelihoods_tol(current_epsilon)
  return(cess(current_log_weights, log_likelihoods) - abc_desired_cess)
}


epsilon_bisection <- function(current_epsilon,
                              all_auxiliary_variables,
                              abc_desired_cess,
                              current_log_weights,
                              all_abc_estimate_log_likelihoods_tol)
{
  current_bisect_epsilon = current_epsilon
  bisect_size = current_bisect_epsilon/2
  direction = -1

  # old_score = cess_score(current_bisect_epsilon,
  #                        all_auxiliary_variables,
  #                        abc_desired_cess,
  #                        current_log_weights,
  #                        all_abc_estimate_log_likelihoods_tol)

  for (i in 1:100)
  {
    new_bisect_epsilon = current_bisect_epsilon + direction*bisect_size

    bisect_size = bisect_size/2

    the_score = cess_score(new_bisect_epsilon,
                           all_auxiliary_variables,
                           abc_desired_cess,
                           current_log_weights,
                           all_abc_estimate_log_likelihoods_tol)

    direction = -sign(the_score)

    current_bisect_epsilon = new_bisect_epsilon

    if (the_score==0)
    {
      break
    }

  }

  return(current_bisect_epsilon)
}


epsilon_doubling <- function(all_auxiliary_variables,
                             abc_desired_cess,
                             current_log_weights,
                             all_abc_estimate_log_likelihoods_tol)
{
  current_epsilon = 1
  current_bisect_epsilon = current_epsilon

  old_score = cess_score(current_bisect_epsilon,
                         all_auxiliary_variables,
                         abc_desired_cess,
                         current_log_weights,
                         all_abc_estimate_log_likelihoods_tol)

  if (old_score>=0)
    return(current_bisect_epsilon)

  for (i in 1:200)
  {
    new_bisect_epsilon = 2*current_bisect_epsilon

    the_score = cess_score(new_bisect_epsilon,
                           all_auxiliary_variables,
                           abc_desired_cess,
                           current_log_weights,
                           all_abc_estimate_log_likelihoods_tol)

    current_bisect_epsilon = new_bisect_epsilon

    if (the_score>=0)
      return(current_bisect_epsilon)

  }

  stop("Doubling method cannot increase abc_tolerance sufficiently to reach desired CESS.")
}
