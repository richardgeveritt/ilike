setup_likelihood_estimator_analytic = function(analytic_evaluate_log_likelihood)
{
  return()
}

make_log_likelihood_estimator = function(model, algorithm)
{
  if (algorithm$likelihood_method=="analytic")
  {
    evaluate_log_likelihood = model$evaluate_log_likelihood
    simulate_auxiliary_variables = function(inputs,data){return(0)}
    setup_likelihood_estimator = function(all_points, all_auxiliary_variables){return(function(point,data,auxiliary_variables){return(evaluate_log_likelihood(point,data))})}
    likelihood_estimator = list(simulate_auxiliary_variables = simulate_auxiliary_variables,
                                setup_likelihood_estimator = setup_likelihood_estimator)

    # virtual double estimate_log_likelihood(const NumericVector &inputs, const NumericVector &data, const List &auxiliary_variables) const=0;
    # virtual void setup_likelihood_estimator(const std::vector<List> &all_auxiliary_variables)=0;
    #
    # estimate_log_likelihood = function(point,data,auxiliary_variables){return(analytic_evaluate_log_likelihood(point,data))}
  }

  return(likelihood_estimator)
}
