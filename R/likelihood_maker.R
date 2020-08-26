
make_likelihood_estimator = function(model, algorithm)
{
  if (algorithm$likelihood_method=="analytic")
  {
    evaluate_log_likelihood = model$evaluate_log_likelihood

    data = model$data

    simulate_auxiliary_variables = function(inputs){return(0)}

    setup_likelihood_estimator = function(all_points, all_auxiliary_variables){return(function(point,auxiliary_variables){return(evaluate_log_likelihood(point,data))})}

    likelihood_estimator = list(simulate_auxiliary_variables = simulate_auxiliary_variables,
                                setup_likelihood_estimator = setup_likelihood_estimator)
  }
  else if (algorithm$likelihood_method=="abc")
  {
    simulate_model = mode$simulate_model

    number_of_simulations = algorithm$number_of_simulations

    evaluate_log_abc_kernel = algorithm$evaluate_log_abc_kernel

    summary_statistic = algorithm$summary_statistic

    abc_tolerances = algorithm$abc_tolerances

    summary_statistic_scaling = algorithm$summary_statistic_scaling

    data = model$data

    summary_data = summary_statistic(data)

    simulate_auxiliary_variables = function(inputs){ abc_simulate_auxiliary_variables(inputs,
                                                                                      number_of_simulations,
                                                                                      simulate_model) }

    setup_likelihood_estimator = function(all_points, all_auxiliary_variables){return(function(point,auxiliary_variables){return(abc_setup_likelihood_estimatorn(all_points,
                                               all_auxiliary_variables,
                                               number_of_simulations,
                                               evaluate_log_abc_kernel,
                                               summary_statistic,
                                               summary_data,
                                               abc_tolerances))})}

    likelihood_estimator = list(simulate_auxiliary_variables = simulate_auxiliary_variables,
                                setup_likelihood_estimator = setup_likelihood_estimator)

  }

  return(likelihood_estimator)
}
