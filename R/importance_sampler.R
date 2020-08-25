

check_model_dimension_consistent_with_simulation = function(stored_parameter_dimension,
                                                            sample_parameter)
{
  sample_parameter_dimension = length(sample_parameter)
  if (is.null(stored_parameter_dimension)) {
    stored_parameter_dimension = sample_parameter_dimension
  }
  else {
    if (stored_parameter_dimension != sample_parameter_dimension) {
      stop("simulated parameter is inconsistent with information in the specified model about the parameter dimension.")
    }
  }

  return(stored_parameter_dimension)
}

check_inputs = function(model,
                        algorithm,
                        is_cpp = FALSE)
{
  # Check that the "inputs" and parameter index are consistent.
  if ( (is.null(model$inputs) ) && (!is.null(model$parameter_index) ) ) {
    # If inputs are not set, simulate them from the proposal.
    if (is_cpp)
    {
      model$inputs = simulate_distribution_cpp(algorithm$simulate_proposal)
    }
    else
    {
      model$inputs = algorithm$simulate_proposal()
    }
  } else if ( (!is.null(model$inputs) ) && (is.null(model$parameter_index) ) ) {
    # If parameter index is not set, simply index all of the inputs.
    model$parameter_index = 1:length(model$inputs)
  } else if ( (!is.null(model$inputs) ) && (!is.null(model$parameter_index) ) ) {
    if (length(model$parameter_index)>length(model$inputs) ) {
      stop("model$parameter_index should be <= the length of model$inputs.")
    }
  } else {
    if (is_cpp)
    {
      model$inputs = simulate_distribution_cpp(algorithm$simulate_proposal)
    }
    else
    {
      model$inputs = algorithm$simulate_proposal()
    }
    model$parameter_index = 1:length(model$inputs)
  }

  if (is_cpp)
  {
    model$parameter_index = model$parameter_index - 1
  }

  return(model)
}


check_IS = function(model,
                    algorithm,
                    is_cpp = FALSE,
                    messages = FALSE)
{
  # Check that any information in the model about the dimension of the parameter is consistent.
  if ( (is.null(model$parameter_dimension) ) && (!is.null(model$parameter_index) ) ) {
    stored_parameter_dimension = length(model$parameter_index)
  } else if ( (!is.null(model$parameter_dimension) ) && (is.null(model$parameter_index) ) ) {
    stored_parameter_dimension = model$parameter_dimension
  } else if ( (!is.null(model$parameter_dimension) ) && (!is.null(model$parameter_index) ) ) {
    if (length(model$parameter_index)==model$parameter_dimension) {
      stored_parameter_dimension = model$parameter_dimension
    } else {
      stop("model$parameter_dimension is not the same as the length of model$parameter_index.")
    }
  } else {
    if (!is.null(model$inputs) )
      stored_parameter_dimension = length(model$inputs)
    else
      stored_parameter_dimension = NULL
  }

  # Simulation cases:
  # Neither prior nor proposal simulator defined.
  # Prior simulator defined; proposal not.
  # Proposal simulator defined; prior not.
  # Prior and proposal simulator defined.
  if (is.null(algorithm$simulate_proposal))
  {
    if (is.null(model$simulate_prior))
    {
      # Neither prior nor proposal simulator defined.
      stop("To use an importance sampler, either model$simulate_prior or algorithm$simulate_proposal need to be specified.")
    }
    else
    {
      # Prior simulator defined; proposal not.
      # In this case, the IS will use the prior as the proposal.
      algorithm$prior_is_proposal = TRUE
      algorithm$simulate_proposal = model$simulate_prior
      if (messages == TRUE)
        print("No method specified for simulating from the proposal. Using the prior as the proposal.")

    }
  }
  else
  {
    if (!is.null(model$simulate_prior))
    {
      # Both prior and proposal simulators are defined, so we should check that their dimensions match.
      if (is_cpp)
      {
        sample_parameter = simulate_distribution_cpp(model$simulate_prior)
      }
      else
      {
        sample_parameter = model$simulate_prior()
      }

      stored_parameter_dimension = check_model_dimension_consistent_with_simulation(stored_parameter_dimension,sample_parameter)
    }

    # We have a standard importance sampler with a different proposal and prior.
    algorithm$prior_is_proposal = FALSE
    if (messaages == TRUE)
      print("Importance points will be simuulate from the specified proposal.")

  }

  # Check all parameter dimensions match.
  if (is_cpp)
  {
    sample_parameter = simulate_distribution_cpp(algorithm$simulate_proposal)
  }
  else
  {
    sample_parameter = algorithm$simulate_proposal()
  }

  stored_parameter_dimension = check_model_dimension_consistent_with_simulation(stored_parameter_dimension,sample_parameter)

  model$parameter_dimension = stored_parameter_dimension

  model = check_inputs(model, algorithm, is_cpp)

  # Check that we can evaluate the prior and the proposal, and that they take as their argument something of dimension length(inputs).
  if (algorithm$prior_is_proposal==FALSE)
  {
    if (is.null(algorithm$evaluate_log_proposal))
    {
      stop("If the prior is not used as the proposal, you need to specify algorithm$evaluate_log_proposal.")
    }
    else
    {
      tryCatch(
        if (is_cpp) {evaluate_log_distribution_cpp(algorithm$evaluate_log_proposal,sample_parameter)} else {algorithm$evaluate_log_proposal(sample_parameter)},
        error = function(e) {stop("algorithm$evaluate_log_proposal generates an error when used on a vector of the dimension of the parameter.")})
    }

    if (is.null(model$evaluate_log_prior))
    {
      stop("If the prior is not used as the proposal, you need to specify algorithm$evaluate_log_proposal.")
    }
    else
    {
      tryCatch(
        if (is_cpp) {evaluate_log_distribution_cpp(algorithm$evaluate_log_prior,sample_parameter)} else {algorithm$evaluate_log_prior(sample_parameter)},
        error = function(e) {stop("model$evaluate_log_prior generates an error when used on a vector of the dimension of the parameter.")})
    }
  }

  # Check that the method for evaluating the likelihood makes sense, and that it takes something of the right dimension.

  if (is.null(model$data))
  {
    stop("model$data not specified.")
  }

  # Check which methods for evaluating the likelihood are available given the specification.
  # Method 1: analytic likelihood.
  # Method 2: standard ABC
  likelihood_methods = matrix(0,2)

  if (!is.null(model$evaluate_log_likelihood))
  {
    likelihood_methods[1] = 1
  }

  if (!is.null(model$simulate_model))
  {
    likelihood_methods[2] = 1
  }

  # Automatically find a method by checking likelihood_methods.
  if (is.null(algorithm$likelihood_method))
  {
    if (likelihood_methods[1] == 1) {
      if (messages == TRUE)
        print("model$evaluate_log_likelihood is specified, so we will use this to evaluate the likelihood.")
      algorithm$likelihood_method = "analytic"
    } else if (likelihood_methods[2] == 1) {
      if (messages == TRUE)
        print("model$simulate_model is specified, so we will use ABC as an approximate likelihood.")
      algorithm$likelihood_method = "abc"
    }

  }


  if (!is.null(algorithm$likelihood_method))
  {
    if (algorithm$likelihood_method=="analytic")
    {
      if (likelihood_methods[1] == 1)
      {
        tryCatch(
          if (is_cpp) {evaluate_log_likelihood_cpp(model$evaluate_log_likelihood,model$inputs,model$data)} else {model$evaluate_log_likelihood(model$inputs,model$data)},
          error = function(e) {stop("model$evaluate_log_likelihood generates an error when used on a vector of dimension model$inputs.")})

        likelihood_estimator = list()
      }
      else
      {
        stop("algorithm$likelihood_method is analytic, but model$evaluate_log_likelihood is not defined.")
      }
    }

    if (algorithm$likelihood_method=="abc")
    {
      if (likelihood_methods[2] == 1)
      {
        tryCatch(
          if (is_cpp) {simulate_model_cpp(model$simulate_model,model$inputs,model$data)} else {model$simulate_model(model$inputs)},
          error = function(e) {stop("model$simulate generates an error when used on a vector of dimension model$inputs.")})

        if (is.null(model$likelihood_options))
        {
          # Need to specify function set up likelihood estimate. Info needs to be stored in a class-like object.
          likelihood_estimator = list(simulate_auxiliary_variables = model$simulate,
                                   setup_likelihood_estimator = ABC$setup_likelihood_estimator)
        }
        else
        {
          stop("Not written yet. Should be checking for: M; epsilon; distance.")
        }

        # Estimator itself needs to take aux variables as an arg, then do the calculation using a method specified in a different file, all of which should be included at the top of this file. Specified after we have sampled all of the thetas, to set thresholds, etc.

      }
      else
      {
        stop("algorithm$likelihood_method is abc, but model$simulate is not defined.")
      }

      # Also should check ABC info.
    }

  }
  else
  {
    stop("No method for evaluating or estimating the likelihood is specified.")
  }

  algorithm$likelihood_estimator = likelihood_estimator

  if (algorithm$likelihood_method!="analytic")
  {
    # Test the generation of the auxiliary variables and store their dimension.
    auxiliary_variables = tryCatch(
      if (is_cpp) {simulate_auxiliary_variables_cpp(algorithm$likelihood_estimator$simulate_auxiliary_variables,model$inputs,model$data)} else {algorithm$likelihood_estimator$simulate_auxiliary_variables(model$inputs,model$data)},
      error = function(e) {stop("algorithm$likelihood_estimator$simulate_auxiliary_variables generates an error when used on a vector of dimension model$inputs.")})

    algorithm$auxiliary_variables_dimension = length(unlist(auxiliary_variables))
  }
  else
  {
    algorithm$auxiliary_variables_dimension = 0
  }


  # Need to check that algorithm$simulate_proposal(2) gives something of the right dimension.

  # Need to check that all of the algorithm$ parts exist.

  list(model=model,algorithm=algorithm)

}


simulate_batch = function(batch_number,
                          num_batch_points,
                          model,
                          algorithm)
{
  # not currently used - needs editing
  proposed_points = algorithm$simulate_proposal(num_batch_points)
  proposed_inputs = t(matrix(rep(model$inputs,num_batch_points),length(model$inputs),num_batch_points))
  proposed_inputs[,model$parameter_index] = proposed_points
  proposed_inputs = lapply(1:num_batch_points,function(i){proposed_inputs[i,]})
  if (!is.null(algorithm$likelihood_estimator))
  {
    proposed_auxiliary_variables = future.apply::future_lapply(proposed_inputs,FUN=function(Input){algorithm$likelihood_estimator$simulate_auxiliary_variables(model$inputs,model$data)},future.seed = TRUE)
  }
}

#' Importance sampler.
#'
#' @param number_of_points The number of importance points.
#' @param model A model.
#' @param algorithm algorithm details.
#' @param max_vector_size If any vector/list has more than this many entries, we need to split the algorithm into multiple batches.
#' @param messages Set to TRUE to have messages printed to console.
#' @return Importance points.
#' @export
importance_sample = function(number_of_points,
                             model,
                             algorithm = list(),
                             max_vector_size = 4e+9,
                             messages = FALSE)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  output = check_IS(model,algorithm)
  model = output$model
  algorithm = output$algorithm

  # Work out the number of batches to simulate.
  # We will store the points and the auxiliary variables in files to avoid memory problems.
  total_dimension = length(model$inputs) + algorithm$auxiliary_variables_dimension
  number_of_batches_minus_one = floor((number_of_points*total_dimension-1)/max_vector_size)

  # simulate batch of parameters and associated aux variables.
  # Save all in file if multiple batches
  # Assume proposal is easy to simulate, and use user-provided function for simulating n times, rather than parallelising.
  # Assume it is expensive to generate the auxiliary variables, and parallelise over this.
  if (number_of_batches_minus_one==0)
  {
    # Everything can be done in one batch, so things are more straightforward.

    # Do the simulation.
    proposed_points = sapply(1:number_of_points,FUN=function(i) {algorithm$simulate_proposal()})
    proposed_inputs = t(matrix(rep(model$inputs,number_of_points),length(model$inputs),number_of_points))
    proposed_inputs[,model$parameter_index] = proposed_points
    proposed_inputs = lapply(1:number_of_points,function(i){proposed_inputs[i,]})

    likelihood_estimator = make_log_likelihood_estimator(model, algorithm)

    proposed_auxiliary_variables = future.apply::future_lapply(proposed_inputs,
                                                               FUN=function(input){likelihood_estimator$simulate_auxiliary_variables(input,model$data)},
                                                               future.seed = TRUE)

    # Now configure the likelihood estimator using all of the simulations, if needed.
    likelihood_estimator$estimate_log_likelihood = tryCatch(likelihood_estimator$setup_likelihood_estimator(proposed_points,proposed_auxiliary_variables),error = function(e) {stop("likelihood_estimator$setup_likelihood_estimator throws an error when used on the proposed points.")})

    log_likelihoods = unlist(future.apply::future_lapply(proposed_inputs,function(i){ likelihood_estimator$estimate_log_likelihood(i, model$data, proposed_auxiliary_variables[[i]]) }))
    # Calculate weights.
    if (algorithm$prior_is_proposal==TRUE)
    {
      log_weights = log_likelihoods
    }
    else
    {
      log_weights = unlist(future.apply::future_lapply(proposed_inputs,function(p){ p = i[model$parameter_index]; model$evaluate_log_prior(p) - algorithm$evaluate_log_proposal(p) })) + log_likelihoods
    }

    Results = list(proposed_points = proposed_points,
                   proposed_auxiliary_variables = proposed_auxiliary_variables,
                   log_weights = log_weights,
                   log_normalising_constant = log_sum_exp(log_weights))

  }
  else
  {
    stop("Importance_sampler not yet set upt to use multiple batches.")

    # We need to use multiple batches since we don't want to store big vectors in memory.
    most_batch_sizes = floor(number_of_points/number_of_batches_minus_one)
    last_batch_size = ((number_of_points*total_dimension)%%max_vector_size)/total_dimension

    # This will store the results to a file.
    lapply(1:number_of_batches_minus_one,FUN=function(i) {simulate_batch(i,most_batch_sizes,model,algorithm); return(NULL)})
    simulate_batch(number_of_batches_minus_one+1,last_batch_size,model,algorithm)
  }

}


#' Importance sampler using cpp functions.
#'
#' @param number_of_points The number of importance points.
#' @param model A model.
#' @param algorithm Algorithm details.
#' @param max_vector_size If any vector/list has more than this many entries, we need to split the algorithm into multiple batches.
#' @param messages Set to TRUE to have messages printed to console.
#' @return Importance points.
#' @export
importance_sample_cpp = function(number_of_points,
                                 model,
                                 algorithm = list(),
                                 max_vector_size = 4e+9,
                                 messages = FALSE)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  output = check_IS(model, algorithm, is_cpp = TRUE, messages = messages)

  model = output$model
  algorithm = output$algorithm

  return(do_importance_sampler_cpp(number_of_points,
                                 model,
                                 algorithm,
                                 max_vector_size))

}
