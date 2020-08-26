evaluate_abc_likelihood = function(simulated_data,
                                   evaluate_log_abc_kernel,
                                   summary_statistic,
                                   summary_data,
                                   abc_tolerances,
                                   summary_statistic_scaling)
{
  return(evaluate_log_abc_kernel(summary_statistic_scaling*summary_statistic(simulated_data),
                                 summary_statistic_scaling*summary_data,
                                 abc_tolerances))
}

abc_simulate_auxiliary_variables = function(inputs,
                                            number_of_simulations,
                                            simulator)
{
  simulations = matrix(0,number_of_simulations,length(inputs))

  repeated_inputs = lapply(1:number_of_simulations,FUN = function(i){inputs})

  simulations = future.apply::future_lapply(1:number_of_simulations,
                              FUN = function(i){simulator(i)},
                              future.seed = TRUE)

  return(simulations)
}

abc_estimate_log_likelihood = function(auxiliary_variables,
                                       evaluate_log_abc_kernel,
                                       summary_statistic,
                                       summary_data,
                                       abc_tolerances,
                                       summary_statistic_scaling)
{
  aux_variables = auxiliary_variables[[1]]

  abc_evaluations = unlist(future.apply::future_lapply(auxiliary_variables,
                                                       FUN = function(i){evaluate_abc_likelihood(i,
                                                                                                 evaluate_log_abc_kernel,
                                                                                                 summary_statistic,
                                                                                                 summary_data,
                                                                                                 abc_tolerances,
                                                                                                 summary_statistic_scaling)},
                                                       future.seed = TRUE))

  return(log_sum_exp(abc_evaluations) - log(length(abc_evaluations)))
}

abc_setup_likelihood_estimator = function(all_points,
                                          all_auxiliary_variables,
                                          number_of_simulations,
                                          evaluate_log_abc_kernel,
                                          summary_statistic,
                                          summary_data,
                                          abc_tolerances)
{
  # Resulting function is the standard ABC estimator: it will take some simulated data variables and find the average of an ABC kernel at these values.

  sum_stat_length = length(summary_data)

  if (length(all_auxiliary_variables)==0)
  {
    summary_statistic_scaling = matrix(1,sum_stat_length)
  }
  else
  {

    data_length = ncol(all_auxiliary_variables[[1]])

    # Use all_auxiliary_variables to find a good scaling for the summary statistics.
    all_sumstats = matrix(0,
                          length(all_auxiliary_variables)*number_of_simulations,
                          sum_stat_length);

    simulated_data_matrix = t(matrix(unlist(all_auxiliary_variables),data_length,length(all_auxiliary_variables)))

    simulated_data_rows = lapply(1:num_points,function(i){simulated_data_matrix[i,]})

    sumstats = future.apply::future_lapply(simulated_data_rows,
                                           FUN=summary_statistic)

    sumstats_vec = sapply(sumstats,c)

    summary_statistic_scaling = 1/apply(sumstats_vec,1,sd,na.rm=TRUE)

  }

  return(function(point,auxiliary_variables){return(abc_estimate_log_likelihood(auxiliary_variables,
                                                                                evaluate_log_abc_kernel,
                                                                                summary_statistics,
                                                                                summary_data,
                                                                                abc_tolerances,
                                                                                summary_statistic_scaling))})

}
