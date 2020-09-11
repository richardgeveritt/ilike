check_parameter_dimension_and_parameter_index = function(model)
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

  return(stored_parameter_dimension)
}


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


check_is_simulate_proposal = function(model, algorithm, stored_parameter_dimension, is_cpp, messages)
{
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

  return(list(model,algorithm))
}


check_is_evaluate_proposal_and_prior = function(model, algorithm, is_cpp)
{
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
        if (is_cpp) {evaluate_log_distribution_cpp(algorithm$evaluate_log_proposal,model$inputs)} else {algorithm$evaluate_log_proposal(model$inputs)},
        error = function(e) {stop("algorithm$evaluate_log_proposal generates an error when used on a vector of the dimension of the parameter.")})
    }

    if (is.null(model$evaluate_log_prior))
    {
      stop("If the prior is not used as the proposal, you need to specify algorithm$evaluate_log_proposal.")
    }
    else
    {
      tryCatch(
        if (is_cpp) {evaluate_log_distribution_cpp(model$evaluate_log_prior,model$inputs)} else {algorithm$evaluate_log_prior(model$inputs)},
        error = function(e) {stop("model$evaluate_log_prior generates an error when used on a vector of the dimension of the parameter.")})
    }
  }
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


check_likelihood_method = function(model, algorithm, is_cpp, messages)
{

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

  # Check that all of the necessary things are specified and working.
  # Provide sensible defaults for things that are not specified.

  if (!is.null(algorithm$likelihood_method))
  {
    ##########################################################
    # Exact likelihood                                       #
    ##########################################################
    # Check:                                                 #
    # - evaluate_log_likelihood                              #
    ##########################################################

    if (algorithm$likelihood_method=="analytic")
    {
      if (likelihood_methods[1] == 1)
      {
        # Check evaluate_log_likelihood.
        tryCatch(
          if (is_cpp) {evaluate_log_likelihood_cpp(model$evaluate_log_likelihood,model$inputs,model$data)} else {model$evaluate_log_likelihood(model$inputs,model$data)},
          error = function(e) {stop("model$evaluate_log_likelihood generates an error when used on a vector of dimension model$inputs.")})
      }
      else
      {
        stop("algorithm$likelihood_method is analytic, but model$evaluate_log_likelihood is not defined.")
      }
    }

    ##########################################################
    # ABC                                                    #
    ##########################################################
    # Check:                                                 #
    # - simulate_model                                       #
    # - number_of_likelihood_particles                       #
    # - get_data_from_simulation                             #
    # - summary_statistics                                  #
    # - summary_statistics_scaling                           #
    # - adapt_summary_statistics_scaling                     #
    # - evaluate_log_abc_kernel                              #
    # - abc_tolerance or abc_desired_cess                    #
    ##########################################################

    if (algorithm$likelihood_method=="abc")
    {
      if (likelihood_methods[2] == 1)
      {

        # Check simulate_model.
        tryCatch(
          if (is_cpp) {simulated = simulate_model_cpp(model$simulate_model, model$inputs, model$data)} else {simulated = model$simulate_model(model$inputs, model$data)},
          error = function(e) {stop("model$simulate generates an error when used on a vector of dimension model$inputs, together with data model$data.")})


        # Check number_of_likelihood_particles.
        if (is.null(algorithm$number_of_likelihood_particles))
        {
          if (messages == TRUE)
            print("algorithm$number_of_likelihood_particles is not set for ABC. Setting it to 1.")
          algorithm$number_of_likelihood_particles = 2
        }


        # Check method for extracting data from simulation.
        if (is.null(algorithm$get_data_from_simulation))
        {
          if (is_cpp)
          {
            if (messages == TRUE)
              print("algorithm$get_data_from_simulation is not set for ABC. Setting it to get the first element of the simulated list.")
            algorithm$get_data_from_simulation = store_get_first_element_of_list_as_numeric_matrix()
          }
          else
          {
            if (is.list(simulated))
            {
              if (messages == TRUE)
                print("algorithm$get_data_from_simulation is not set for ABC. Setting it to get the first element of the simulated list.")
              algorithm$get_data_from_simulation = function(s){return(s[[1]])}
            }
            else
            {
              algorithm$get_data_from_simulation = function(s){return(s)}
            }

          }
        }
        tryCatch(
          if (is_cpp) {simulated_data = get_data_from_simulation_cpp(algorithm$get_data_from_simulation, simulated)} else {simulated_data = algorithm$get_data_from_simulation(simulated)},
          error = function(e) {stop("algorithm$get_data_from_simulation generates an error when used on a simulation.")})


        # Check summary_statistics.
        if (is.null(algorithm$summary_statistics))
        {
          if (messages == TRUE)
            print("algorithm$summary_statistics is not set for ABC. Setting it to be the identity.")
          if (is_cpp)
          {
            algorithm$summary_statistics = store_identity_statistic()
          }
          else
          {
            algorithm$summary_statistics = function(d){return(d)}
          }
        }
        tryCatch(
          if (is_cpp) {summary_simulated = summary_statistics_cpp(algorithm$summary_statistics, simulated_data)} else {summary_simulated = algorithm$summary_statistics(simulated_data)},
          error = function(e) {stop("algorithm$summary_statistics generates an error when used on simulated data.")})
        tryCatch(
          if (is_cpp) {summary_observed = summary_statistics_cpp(algorithm$summary_statistics, model$data)} else {summary_observed = algorithm$summary_statistics(model$data)},
          error = function(e) {stop("algorithm$summary_statistics generates an error when used on observed data.")})


        # Check summary_statistics_scaling.
        if (is.null(algorithm$summary_statistics_scaling))
        {
          if (messages == TRUE)
            print("algorithm$summary_statistics_scaling is not set for ABC. Setting it to ones.")
          algorithm$summary_statistics_scaling = matrix(1,length(summary_simulated))
        }
        if (length(algorithm$summary_statistics_scaling)!=length(summary_simulated))
          stop("algorithm$summary_statistics_scaling is not the same length as the summary statistic vector.")


        # Check adapt_summary_statistics_scaling.
        if (is.null(algorithm$adapt_summary_statistics_scaling))
        {
          if (messages == TRUE)
            print("algorithm$adapt_summary_statistics_scaling is not set for ABC. Setting it to be TRUE.")
          algorithm$adapt_summary_statistics_scaling = TRUE
        }


        # Check evaluate_log_abc_kernel.
        if (is.null(algorithm$evaluate_log_abc_kernel))
        {
          if (is.null(algorithm$abc_kernel_method))
          {
            if (messages == TRUE)
              print("algorithm$evaluate_log_abc_kernel is not set for ABC. Setting it to uniform, with L2 norm.")
            algorithm$abc_kernel_method = "L2_uniform"
          }

          if (algorithm$abc_kernel_method=="L1_uniform")
          {
            if (is_cpp)
            {
              algorithm$evaluate_log_abc_kernel = store_L1_uniform_evaluate_log_abc_kernel()
            }
            else
            {
              algorithm$evaluate_log_abc_kernel = function(s_s,s_o,tol){return(Lp_uniform_evaluate_log_abc_kernel(s_s,s_o,tol,1))}
            }
          }
          else if (algorithm$abc_kernel_method=="L2_uniform")
          {
            if (is_cpp)
            {
              algorithm$evaluate_log_abc_kernel = store_L2_uniform_evaluate_log_abc_kernel()
            }
            else
            {
              algorithm$evaluate_log_abc_kernel = function(s_s,s_o,tol){return(Lp_uniform_evaluate_log_abc_kernel(s_s,s_o,tol,2))}
            }
          }
          else if (algorithm$abc_kernel_method=="Linf_uniform")
          {
            if (is_cpp)
            {
              algorithm$evaluate_log_abc_kernel = store_Linf_uniform_evaluate_log_abc_kernel()
            }
            else
            {
              algorithm$evaluate_log_abc_kernel = function(s_s,s_o,tol){return(Linf_uniform_evaluate_log_abc_kernel(s_s,s_o,tol))}
            }
          }
          else if (algorithm$abc_kernel_method=="gaussian")
          {
            if (is_cpp)
            {
              algorithm$evaluate_log_abc_kernel = store_gaussian_evaluate_log_abc_kernel()
            }
            else
            {
              algorithm$evaluate_log_abc_kernel = function(s_s,s_o,tol){return(gaussian_evaluate_log_abc_kernel(s_s,s_o,tol))}
            }
          }
          else
          {
            stop("algorithm$abc_kernel_method set to an invalid option.")
          }
        }
        else
        {

          if (messages == TRUE)
            print("algorithm$evaluate_log_abc_kernel is set, so ignoring algorithm$abc_kernel_method.")

          tryCatch(
            if (is_cpp) {evaluate_log_abc_kernel_cpp(algorithm$evaluate_log_abc_kernel, algorithm$summary_statistics_scaling*summary_simulated, algorithm$summary_statistics_scaling*summary_simulated, 1)} else {algorithm$evaluate_log_abc_kernel(algorithm$summary_statistics_scaling*summary_simulated, algorithm$summary_statistics_scaling*summary_observed, 1)},
            error = function(e) {stop("algorithm$evaluate_log_abc_kernel generates an error when used on simulated data and model$data.")})
        }


        # Check abc_tolerance - need either this or a percentage of effectively independent simulations to keep.
        # Need also to add this into setup llhd.
        # If not set, set to keep 100 effective points.
        # There are three options:
        # 1. Set a fixed ABC tolerance.
        # 2. Set a desired CESS.
        # 3. Estimate it from the data.
        # Only one of these should be set. Throw an error if more than one of them is.
        # If it is estimated:
        # - the other options are not possible
        # - we need to have some inputs index set. If no inputs are set, it should be tagged onto the end, and some default prior set
        if (is.null(algorithm$abc_tolerance))
        {
          algorithm$adapt_abc_tolerance_to_cess = TRUE
          #algorithm$abc_tolerance = -1

          if (is.null(algorithm$abc_desired_cess))
          {
            if (messages == TRUE)
              print("algorithm$abc_desired_cess and algorithm$abc_tolerance are not set for ABC. Setting algorithm$abc_desired_cess to 100.")
            algorithm$abc_desired_cess = min(algorithm$number_of_points, 100)
          }
          else
          {
            if (messages == TRUE)
              print("algorithm$abc_tolerance not set for ABC. Using algorithm$abc_desired_cess instead.")

            if (algorithm$abc_desired_cess<1)
            {
              if (messages == TRUE)
                print("algorithm$abc_desired_cess below 1. Set to 1.")
              algorithm$abc_desired_cess = 1

            }
            else if (algorithm$abc_desired_cess>algorithm$number_of_points)
            {
              if (messages == TRUE)
                print("algorithm$abc_desired_cess above algorithm$number_of_points. Set to algorithm$number_of_points.")
              algorithm$abc_desired_cess = algorithm$number_of_points
            }

          }
        }
        else
        {
          algorithm$adapt_abc_tolerance_to_cess = FALSE

          if (is.null(algorithm$abc_desired_cess))
          {
            if (messages == TRUE)
              print("algorithm$abc_desired_cess not set for ABC. Using algorithm$abc_tolerance instead.")

            if (length(algorithm$abc_tolerance)!=1)
              stop("algorithm$abc_tolerance is not a scalar.")

            if (algorithm$abc_tolerance<0)
              stop("algorithm$abc_tolerance is negative.")
          }

          stop("algorithm$abc_tolerance and algorithm$abc_desired_cess both defined. This is ambiguous, and you need to get rid of one of them.")
        }

      }
      else
      {
        stop("algorithm$likelihood_method is abc, but model$simulate is not defined.")
      }

    }

  }
  else
  {
    stop("No method for evaluating or estimating the likelihood is specified.")
  }

  # if (algorithm$likelihood_method!="analytic")
  # {
  #   # Test the generation of the auxiliary variables and store their dimension.
  #   auxiliary_variables = tryCatch(
  #     if (is_cpp) {simulate_auxiliary_variables_cpp(algorithm$likelihood_estimator$simulate_auxiliary_variables,model$inputs,model$data)} else {algorithm$likelihood_estimator$simulate_auxiliary_variables(model$inputs,model$data)},
  #     error = function(e) {stop("algorithm$likelihood_estimator$simulate_auxiliary_variables generates an error when used on a vector of dimension model$inputs.")})
  #
  #   algorithm$auxiliary_variables_dimension = length(unlist(auxiliary_variables))
  # }
  # else
  # {
  #   algorithm$auxiliary_variables_dimension = 0
  # }
  #
  # algorithm$likelihood_estimator = likelihood_estimator

  return(list(model=model, algorithm=algorithm))
}


check_is = function(model,
                    algorithm,
                    is_cpp = FALSE,
                    messages = FALSE)
{

  # Check data is specified.
  if ( (is.null(algorithm$future_type)) && (is.null(algorithm$internal_future_type)) )
  {
    if (messages == TRUE)
      print("algorithm$future_type is not set. Setting it to sequential.")
    future::plan(list("sequential"))
    likelihood_use_future <<- FALSE
  }
  else if ( (is.null(algorithm$future_type)) && (!is.null(algorithm$internal_future_type)) )
  {
    tryCatch(future::plan(list("sequential", algorithm$internal_future_type)),
             error = function(e) {stop("algorithm$internal_future_type is defined but of an invalid type.")})
    likelihood_use_future <<- TRUE
  }
  else if ( (!is.null(algorithm$future_type)) && (is.null(algorithm$internal_future_type)) )
  {
    tryCatch(future::plan(list(algorithm$future_type)),
             error = function(e) {stop("algorithm$future_type is defined but of an invalid type.")})
    likelihood_use_future <<- FALSE
  }
  else if ( (!is.null(algorithm$future_type)) && (!is.null(algorithm$internal_future_type)) )
  {
    tryCatch(future::plan(list(algorithm$future_type, algorithm$internal_future_type)),
             error = function(e) {stop("algorithm$future_type and algorithm$internal_future_type is defined but at least one of them is an invalid type.")})
    likelihood_use_future <<- TRUE
  }

  # Check data is specified.
  if (is.null(algorithm$number_of_points))
  {
    if (messages == TRUE)
      print("algorithm$number_of_points is not set. Setting it to 10000.")
    algorithm$number_of_points = 10000
  }

  # Check consistency of optional arguments.
  stored_parameter_dimension = check_parameter_dimension_and_parameter_index(model)

  # Check mechanism for simulating proposal and make sure it is consistent with other information about the parameter dimension.
  output = check_is_simulate_proposal(model, algorithm, stored_parameter_dimension, is_cpp, messages)
  model = output[[1]]
  algorithm = output[[2]]

  # Check model$inputs (if specified) is consistent with simulated parameters.
  model = check_inputs(model, algorithm, is_cpp)

  # Check methods for evaluating proposal and prior.
  check_is_evaluate_proposal_and_prior(model, algorithm)

  # Check data is specified.
  if (is.null(model$data))
  {
    stop("model$data not specified.")
  }
  else
  {
    model$data = as.matrix(model$data)
  }

  # Check that the method for evaluating the likelihood makes sense, and that it takes something of the right dimension.
  output = check_likelihood_method(model, algorithm, is_cpp, messages)
  model = output$model
  algorithm = output$algorithm

  # Things to do:
  # Need to check that algorithm$simulate_proposal(2) gives something of the right dimension.
  # Need to check that all of the algorithm$ parts exist.

  list(model=model, algorithm=algorithm)

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
}

#' Importance sampler.
#'
#' @param model A model.
#' @param algorithm algorithm details.
#' @param max_vector_size If any vector/list has more than this many entries, we need to split the algorithm into multiple batches.
#' @param messages Set to TRUE to have messages printed to console.
#' @return Importance points.
#' @export
importance_sample = function(model,
                             algorithm,
                             max_vector_size = 4e+9,
                             messages = FALSE)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  output = check_is(model,algorithm)
  model = output$model
  algorithm = output$algorithm

  number_of_points = algorithm$number_of_points

  # Work out the number of batches to simulate.
  # We will store the points and the auxiliary variables in files to avoid memory problems.
  total_dimension = length(model$inputs)# + algorithm$auxiliary_variables_dimension
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

    likelihood_estimator = make_likelihood_estimator(model, algorithm)

    proposed_auxiliary_variables = future.apply::future_lapply(proposed_inputs,
                                                               FUN=function(input){likelihood_estimator$simulate_auxiliary_variables(input)},
                                                               future.seed = TRUE)

    # Now configure the likelihood estimator using all of the simulations, if needed.
    likelihood_estimator$estimate_log_likelihood = likelihood_estimator$setup_likelihood_estimator(proposed_points,proposed_auxiliary_variables)#tryCatch(likelihood_estimator$setup_likelihood_estimator(proposed_points,proposed_auxiliary_variables),error = function(e) {stop("likelihood_estimator$setup_likelihood_estimator throws an error when used on the proposed points.")})

    log_likelihoods = unlist(future.apply::future_lapply(1:length(proposed_inputs),function(i){ likelihood_estimator$estimate_log_likelihood(proposed_inputs[[i]], proposed_auxiliary_variables[[i]]) }))

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
    stop("Importance_sampler not yet set up to use multiple batches.")

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
#' @param model A model.
#' @param algorithm Algorithm details.
#' @param max_vector_size If any vector/list has more than this many entries, we need to split the algorithm into multiple batches.
#' @param messages Set to TRUE to have messages printed to console.
#' @return Importance points.
#' @export
importance_sample_cpp = function(model,
                                 algorithm,
                                 max_vector_size = 4e+9,
                                 messages = FALSE)
{
  # Check that inputs make sense.
  # Need also to:
  # - get dimensions of inputs/parameters
  # - get dimensions of aux variables
  # - make sure indexing of inputs/parameters is stored if necessary
  # - make sure prior is in standard format and store it
  output = check_is(model, algorithm, is_cpp = TRUE, messages = messages)
  model = output$model
  algorithm = output$algorithm

  return(do_importance_sampler_cpp(model,
                                   algorithm,
                                   max_vector_size))

}
