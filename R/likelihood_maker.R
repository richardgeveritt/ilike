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
    simulate_model = model$simulate_model

    get_data_from_simulation = algorithm$get_data_from_simulation

    number_of_likelihood_particles = algorithm$number_of_likelihood_particles

    evaluate_log_abc_kernel = algorithm$evaluate_log_abc_kernel

    summary_statistics = algorithm$summary_statistics

    abc_tolerance = algorithm$abc_tolerance

    abc_desired_cess = algorithm$abc_desired_cess

    summary_statistics_scaling = algorithm$summary_statistics_scaling

    adapt_abc_tolerance_to_cess = algorithm$adapt_abc_tolerance_to_cess

    adapt_summary_statistics_scaling = algorithm$adapt_summary_statistics_scaling

    data = model$data

    summary_data = summary_statistics(data)

    simulate_auxiliary_variables = function(inputs){ abc_simulate_auxiliary_variables(inputs,
                                                                                      number_of_likelihood_particles,
                                                                                      simulate_model,
                                                                                      data) }

    setup_likelihood_estimator = function(all_points, all_auxiliary_variables){return(abc_is_setup_likelihood_estimator(all_points,
                                              all_auxiliary_variables,
                                              number_of_likelihood_particles,
                                              get_data_from_simulation,
                                              evaluate_log_abc_kernel,
                                              summary_statistics,
                                              summary_data,
                                              abc_tolerance,
                                              abc_desired_cess,
                                              summary_statistics_scaling,
                                              adapt_abc_tolerance_to_cess,
                                              adapt_summary_statistics_scaling))}

    likelihood_estimator = list(simulate_auxiliary_variables = simulate_auxiliary_variables,
                                setup_likelihood_estimator = setup_likelihood_estimator)

  }

  return(likelihood_estimator)
}
