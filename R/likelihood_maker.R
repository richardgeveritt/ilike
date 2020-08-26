setup_likelihood_estimator_analytic = function(analytic_evaluate_log_likelihood)
{
  return()
}

make_likelihood_estimator = function(model, algorithm)
{
  if (algorithm$likelihood_method=="analytic")
  {
    evaluate_log_likelihood = model$evaluate_log_likelihood
    simulate_auxiliary_variables = function(inputs,data){return(0)}
    setup_likelihood_estimator = function(all_points, all_auxiliary_variables){return(function(point,data,auxiliary_variables){return(evaluate_log_likelihood(point,data))})}
    likelihood_estimator = list(simulate_auxiliary_variables = simulate_auxiliary_variables,
                                setup_likelihood_estimator = setup_likelihood_estimator)
  }
  else if (algorithm$likelihood_method=="abc")
  {

  }

  return(likelihood_estimator)
}
