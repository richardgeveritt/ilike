split_string_at_comma_ignoring_parentheses <- function(input_string)
{
  result <- vector()
  current_word <- ""
  parentheses_level <- 0

  for (char in strsplit(input_string, "")[[1]]) {
    if (char == "," && parentheses_level == 0) {
      result <- c(result, current_word)
      current_word <- ""
    } else {
      current_word <- paste0(current_word, char)

      if (char == "(") {
        parentheses_level <- parentheses_level + 1
      } else if (char == ")") {
        parentheses_level <- parentheses_level - 1
      }
    }
  }

  result <- c(result, current_word)

  # Trim leading and trailing whitespace from each split substring
  result <- trimws(result)

  return(result)
}

my_julia_source = function(filename,
                           julia_required_libraries=c())
{
  # Read the content of the Julia file
  file_content <- readLines(filename)

  # Define a regular expression pattern to identify function definitions
  function_pattern <- "^function\\s+(\\w+)"

  # Initialize variables
  current_function <- NULL
  current_lines <- NULL

  output = list()
  list_names = c()

  # Loop through each line in the file
  for (line in file_content) {
    # Check if the line matches the function pattern
    if (grepl(function_pattern, line)) {
      # If a function is already being processed, write it to a file
      if (!is.null(current_function)) {

        function_name <- strsplit(current_function, "\\(")[[1]][1]
        #assign(function_name,JuliaConnectoR::juliaEval(paste(current_lines,collapse="\n")))
        assign(function_name,JuliaCall::julia_eval(paste(current_lines,collapse="\n")))
        output = append(output,eval(parse(text=function_name)))
        list_names = c(list_names,function_name)
        #writeLines(current_lines, paste0(current_function, ".jl"))
      }
      else {
        #JuliaConnectoR::juliaEval(paste(current_lines,collapse="\n"))
        #JuliaCall::julia_eval(paste(current_lines,collapse="\n"))
        if (length(julia_required_libraries)>0)
        {
          for (i in 1:length(julia_required_libraries))
          {
            JuliaCall::julia_install_package(julia_required_libraries[i])
            JuliaCall::julia_library(julia_required_libraries[i])
          }
        }
      }

      # Extract the function name from the line
      current_function <- sub(function_pattern, "\\1", line)
      # Reset the lines buffer
      current_lines <- NULL
    }

    # Add the current line to the lines buffer
    current_lines <- c(current_lines, line)
  }

  # Write the last function to a file
  if (!is.null(current_function)) {

    function_name <- strsplit(current_function, "\\(")[[1]][1]
    #assign(function_name,JuliaConnectoR::juliaEval(paste(current_lines,collapse="\n")))
    assign(function_name,JuliaCall::julia_eval(paste(current_lines,collapse="\n")))
    output = append(output,eval(parse(text=function_name)))
    list_names = c(list_names,function_name)

    #writeLines(current_lines, paste0(current_function, ".jl"))
  }

  names(output) = list_names
  return(output)
}

ilike_parse <- function(input,
                        parameter_list = list())
{
  required_args <- formula.tools::get.vars(parse(text=input))

  an.error.occured <- FALSE
  tryCatch( { required_args <- formula.tools::get.vars(parse(text=input)) }
            , error = function(e) {an.error.occured <<- TRUE})
  if (an.error.occured)
  {
    stop(paste('Error in formula.tools::get.vars function. Maybe you used an equals sign within the function specification (e.g. using log=TRUE). This is not yet supported.',sep=""))
  }

  parameter_arguments = c()

  if (length(required_args)>0)
  {
    for (i in 1:length(required_args))
    {
      h = required_args[i]

      if ( (nchar(h)>1) && (substr(h,1,1)=="p") && (!grepl("\\D",substr(h,2,nchar(h)))) )
      {
        parameter_number = as.numeric(substr(h,2,nchar(h)))
        if (parameter_number<=length(parameter_list))
        {
          parameter_arguments = c(parameter_arguments,h)
          do.call("<-",list(h, parameter_list[[parameter_number]]))
        }
        else
        {
          stop(paste("In call ",input,", parameter ",as.numeric(substr(h,2,nchar(h)))," not found in parameter_list (parameter_list has length ",length(parameter_list),").",sep=""))
        }
      }
    }

    # If p1, etc are in the list, remove from arg list, set p1=, etc. Then should be used by function.
    required_args = setdiff(required_args,parameter_arguments)
  }

  #eval(parse(text = paste('result_function <- function(', paste(required_args,collapse=","), ') { return(' , input , ')}', sep='')))
  result_function_text = paste(' <- function(', paste(required_args,collapse=","), ') { return(' , input , ')}', sep='')
  return(list(result_function_text,required_args))
}

factor_processing = function(factor_number,blocks,block_name,prior_function_types,custom_likelihood_function_types,sbi_likelihood_function_types,linear_gaussian_data_model_types,nonlinear_gaussian_data_model_types,other_likelihood_function_types,line_counter)
{
  # Is this a continuation of the current factor, or a new one?

  all_names = c(prior_function_types,custom_likelihood_function_types,sbi_likelihood_function_types,linear_gaussian_data_model_types,nonlinear_gaussian_data_model_types,other_likelihood_function_types)

  # Get the current factor info.
  if ("factor" %in% names(blocks))
  {
    current_factor_info = blocks[["factor"]][[factor_number]]

    current_factor_names = names(current_factor_info)

    if (block_name %in% prior_function_types)
    {
      if ("prior" %in% names(current_factor_info))
      {
        # factor is complete
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("evaluate_log_prior" %in% names(current_factor_info)) || ("simulate_prior" %in% names(current_factor_info)) || ("evaluate_gradient_log_prior" %in% names(current_factor_info)) || ("evaluate_second_gradient_log_prior" %in% names(current_factor_info)) )
      {
        # factor is complete
        if ( ("evaluate_log_prior" %in% names(current_factor_info)) && (block_name=="evaluate_log_prior") )
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
        }
        else if ( ("simulate_prior" %in% names(current_factor_info)) && (block_name=="simulate_prior") )
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
        }
        else if ( ("evaluate_gradient_log_prior" %in% names(current_factor_info)) && (block_name=="evaluate_gradient_log_prior") )
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
        }
        else if ( ("evaluate_second_gradient_log_prior" %in% names(current_factor_info)) && (block_name=="evaluate_second_gradient_log_prior") )
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
        }
      }

      # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
      all_other_names = setdiff(all_names,prior_function_types)
      for (i in 1:length(current_factor_names))
      {
        if ( current_factor_names[i] %in% all_other_names)
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
          break
        }
      }
    }
    else if (block_name %in% custom_likelihood_function_types)
    {
      # factor is complete
      if ( ("evaluate_log_likelihood" %in% names(current_factor_info)) && (block_name=="evaluate_log_likelihood") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("evaluate_gradient_log_likelihood" %in% names(current_factor_info)) && (block_name=="evaluate_gradient_log_likelihood") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("evaluate_second_gradient_log_likelihood" %in% names(current_factor_info)) && (block_name=="evaluate_second_gradient_log_likelihood") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }

      # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
      all_other_names = setdiff(all_names,custom_likelihood_function_types)
      for (i in 1:length(current_factor_names))
      {
        if ( current_factor_names[i] %in% all_other_names)
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
          break
        }
      }
    }
    # else if (block_name %in% file_likelihood_function_types)
    # {
    #   # factor is complete
    #   if ( ("importance_sample" %in% names(current_factor_info)) && (block_name=="importance_sample") )
    #   {
    #     print_factor_info(factor_number,blocks,line_counter-1)
    #     factor_number = factor_number + 1
    #   }
    #   else if ( ("smc_mcmc_move" %in% names(current_factor_info)) && (block_name=="smc_mcmc_move") )
    #   {
    #     print_factor_info(factor_number,blocks,line_counter-1)
    #     factor_number = factor_number + 1
    #   }
    #
    #   # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
    #   all_other_names = setdiff(all_names,file_likelihood_function_types)
    #   for (i in 1:length(current_factor_names))
    #   {
    #     if ( current_factor_names[i] %in% all_other_names)
    #     {
    #       print_factor_info(factor_number,blocks,line_counter-1)
    #       factor_number = factor_number + 1
    #       break
    #     }
    #   }
    # }
    else if (block_name %in% sbi_likelihood_function_types)
    {
      # factor is complete
      if ( ("simulate_data_model" %in% names(current_factor_info)) && (block_name=="simulate_data_model") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("sbi_likelihood" %in% names(current_factor_info)) && (block_name=="sbi_likelihood") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("summary_statistics" %in% names(current_factor_info)) && (block_name=="summary_statistics") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }

      # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
      all_other_names = setdiff(all_names,sbi_likelihood_function_types)
      for (i in 1:length(current_factor_names))
      {
        if ( current_factor_names[i] %in% all_other_names)
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
          break
        }
      }
    }
    else if (block_name %in% linear_gaussian_data_model_types)
    {
      # factor is complete
      if ( ("linear_gaussian_data_model" %in% names(current_factor_info)) && (block_name=="linear_gaussian_data_model") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("linear_gaussian_data_matrix" %in% names(current_factor_info)) && (block_name=="linear_gaussian_data_matrix") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("linear_gaussian_data_covariance" %in% names(current_factor_info)) && (block_name=="linear_gaussian_data_covariance") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }

      # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
      all_other_names = setdiff(all_names,linear_gaussian_data_model_types)
      for (i in 1:length(current_factor_names))
      {
        if ( current_factor_names[i] %in% all_other_names)
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
          break
        }
      }
    }
    else if (block_name %in% nonlinear_gaussian_data_model_types)
    {
      # factor is complete
      if ( ("nonlinear_gaussian_data_model" %in% names(current_factor_info)) && (block_name=="nonlinear_gaussian_data_model") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("nonlinear_gaussian_data_function" %in% names(current_factor_info)) && (block_name=="nonlinear_gaussian_data_function") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("nonlinear_gaussian_data_covariance" %in% names(current_factor_info)) && (block_name=="nonlinear_gaussian_data_covariance") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }

      # If any of the existing parts of the factor are of a completely different type to the block_name then make a new factor.
      all_other_names = setdiff(all_names,nonlinear_gaussian_data_model_types)
      for (i in 1:length(current_factor_names))
      {
        if ( current_factor_names[i] %in% all_other_names)
        {
          print_factor_info(factor_number,blocks,line_counter-1)
          factor_number = factor_number + 1
          break
        }
      }
    }
    else if (block_name %in% other_likelihood_function_types)
    {
      print_factor_info(factor_number,blocks,line_counter-1)
      factor_number = factor_number + 1
    }
    else
    {
      stop(paste("Invalid block: ",block_name,'.',sep=""))
    }
  }
  else
  {
    factor_number = factor_number + 1
  }

  return(factor_number)
}

transition_model_processing = function(transition_model_number,blocks,block_name,ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_types,line_counter)
{
  # Is this a continuation of the current transition_model, or a new one?

  all_names = c(ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_typess)

  # Get the current transition_model info.
  if ("transition_model" %in% names(blocks))
  {
    current_transition_model_info = blocks[["transition_model"]][[transition_model_number]]

    current_transition_model_names = names(current_transition_model_info)

    if (block_name %in% ilike_transition_model_function_types)
    {
      # transition_model is complete
      if ( ("transition_model" %in% names(current_transition_model_info)) && (block_name=="transition_model") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }

      # If any of the existing parts of the transition_model are of a completely different type to the block_name then make a new transition_model.
      all_other_names = setdiff(all_names,ilike_likelihood_function_types)
      for (i in 1:length(current_transition_model_names))
      {
        if ( current_transition_model_names[i] %in% all_other_names)
        {
          print_transition_model_info(transition_model_number,blocks,line_counter-1)
          transition_model_number = transition_model_number + 1
          break
        }
      }
    }
    else if (block_name %in% custom_transition_model_function_types)
    {
      # transition_model is complete
      if ( ("evaluate_log_transition_model" %in% names(current_transition_model_info)) && (block_name=="evaluate_log_transition_model") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }
      else if ( ("simulate_transition_model" %in% names(current_transition_model_info)) && (block_name=="simulate_transition_model") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }

      # If any of the existing parts of the transition_model are of a completely different type to the block_name then make a new transition_model.
      all_other_names = setdiff(all_names,custom_likelihood_function_types)
      for (i in 1:length(current_transition_model_names))
      {
        if ( current_transition_model_names[i] %in% all_other_names)
        {
          print_transition_model_info(transition_model_number,blocks,line_counter-1)
          transition_model_number = transition_model_number + 1
          break
        }
      }
    }
    else if (block_name %in% linear_gaussian_transition_model_types)
    {
      # transition_model is complete
      if ( ("linear_gaussian_transition_model" %in% names(current_transition_model_info)) && (block_name=="linear_gaussian_transition_model") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }
      else if ( ("linear_gaussian_transition_matrix" %in% names(current_transition_model_info)) && (block_name=="linear_gaussian_transition_matrix") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }
      else if ( ("linear_gaussian_transition_covariance" %in% names(current_transition_model_info)) && (block_name=="linear_gaussian_transition_covariance") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }

      # If any of the existing parts of the transition_model are of a completely different type to the block_name then make a new transition_model.
      all_other_names = setdiff(all_names,linear_gaussian_transition_model_types)
      for (i in 1:length(current_transition_model_names))
      {
        if ( current_transition_model_names[i] %in% all_other_names)
        {
          print_transition_model_info(transition_model_number,blocks,line_counter-1)
          transition_model_number = transition_model_number + 1
          break
        }
      }
    }
    else if (block_name %in% nonlinear_gaussian_transition_model_types)
    {
      # transition_model is complete
      if ( ("nonlinear_gaussian_transition_model" %in% names(current_transition_model_info)) && (block_name=="nonlinear_gaussian_transition_model") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }
      else if ( ("nonlinear_gaussian_transition_function" %in% names(current_transition_model_info)) && (block_name=="nonlinear_gaussian_transition_function") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }
      else if ( ("nonlinear_gaussian_transition_covariance" %in% names(current_transition_model_info)) && (block_name=="nonlinear_gaussian_transition_covariance") )
      {
        print_transition_model_info(transition_model_number,blocks,line_counter-1)
        transition_model_number = transition_model_number + 1
      }

      # If any of the existing parts of the transition_model are of a completely different type to the block_name then make a new transition_model.
      all_other_names = setdiff(all_names,nonlinear_gaussian_transition_model_types)
      for (i in 1:length(current_transition_model_names))
      {
        if ( current_transition_model_names[i] %in% all_other_names)
        {
          print_transition_model_info(transition_model_number,blocks,line_counter-1)
          transition_model_number = transition_model_number + 1
          break
        }
      }
    }
    else if (block_name %in% other_likelihood_function_types)
    {
      print_transition_model_info(transition_model_number,blocks,line_counter-1)
      transition_model_number = transition_model_number + 1
    }
    else
    {
      stop(paste("Invalid block: ",block_name,'.',sep=""))
    }
  }
  else
  {
    transition_model_number = transition_model_number + 1
  }

  return(transition_model_number)
}

potential_function_processing = function(transition_model_number,blocks,block_name,custom_potential_function_types,ilike_potential_function_types,line_counter)
{
  # Is this a continuation of the current potential_function, or a new one?

  all_names = c(custom_potential_function_types,ilike_potential_function_types)

  # Get the current potential_function info.
  if ("potential_function" %in% names(blocks))
  {
    current_potential_function_info = blocks[["potential_function"]][[potential_function_number]]

    current_potential_function_names = names(current_potential_function_info)

    if (block_name %in% ilike_potential_function_function_types)
    {
      # potential_function is complete
      if ( ("potential_function" %in% names(current_potential_function_info)) && (block_name=="potential_function") )
      {
        print_potential_function_info(potential_function_number,blocks,line_counter-1)
        potential_function_number = potential_function_number + 1
      }

      # If any of the existing parts of the potential_function are of a completely different type to the block_name then make a new potential_function.
      all_other_names = setdiff(all_names,ilike_likelihood_function_types)
      for (i in 1:length(current_potential_function_names))
      {
        if ( current_potential_function_names[i] %in% all_other_names)
        {
          print_potential_function_info(potential_function_number,blocks,line_counter-1)
          potential_function_number = potential_function_number + 1
          break
        }
      }
    }
    if (block_name %in% custom_potential_function_function_types)
    {
      # potential_function is complete
      if ( ("evaluate_log_potential_function" %in% names(current_potential_function_info)) && (block_name=="evaluate_log_potential_function") )
      {
        print_potential_function_info(potential_function_number,blocks,line_counter-1)
        potential_function_number = potential_function_number + 1
      }

      # If any of the existing parts of the potential_function are of a completely different type to the block_name then make a new potential_function.
      all_other_names = setdiff(all_names,custom_likelihood_function_types)
      for (i in 1:length(current_potential_function_names))
      {
        if ( current_potential_function_names[i] %in% all_other_names)
        {
          print_potential_function_info(potential_function_number,blocks,line_counter-1)
          potential_function_number = potential_function_number + 1
          break
        }
      }
    }
    else
    {
      stop(paste("Invalid block: ",block_name,'.',sep=""))
    }
  }
  else
  {
    potential_function_number = potential_function_number + 1
  }

  return(potential_function_number)
}

data_processing = function(data_number,line_counter)
{
  if (data_number!=0)
  {
    stop(paste("Invalid file: line ",line_counter,', data already specified.',sep=""))
  }
  data_number = data_number + 1
  return(data_number)
}

importance_proposal_processing = function(importance_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current importance_proposal, or a new one?

  # Get the current importance_proposal info.
  if ("importance_proposal" %in% names(blocks))
  {
    current_importance_proposal_info = blocks[["importance_proposal"]][[importance_proposal_number]]

    if ("importance_proposal" %in% names(current_importance_proposal_info))
    {
      # importance_proposal is complete
      print_importance_proposal_info(importance_proposal_number,blocks,line_counter-1)
      importance_proposal_number = importance_proposal_number + 1
    }
    else if ( ("evaluate_log_importance_proposal" %in% names(current_importance_proposal_info)) || ("simulate_importance_proposal" %in% names(current_importance_proposal_info)) )
    {
      # importance_proposal is complete
      if ( ("evaluate_log_importance_proposal" %in% names(current_importance_proposal_info)) && ("simulate_importance_proposal" %in% names(current_importance_proposal_info)) )
      {
        print_importance_proposal_info(importance_proposal_number,blocks,line_counter-1)
        importance_proposal_number = importance_proposal_number + 1
      }

      if (block_name=="importance_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', importance_proposal specified, but previous importance proposal block is incomplete (did you specify both evaluate_log_importance_proposal and simulate_importance_proposal?)',sep=""))
      }

    }
  }
  else
  {
    importance_proposal_number = importance_proposal_number + 1
  }

  return(importance_proposal_number)
}

mh_proposal_processing = function(mh_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current mh_proposal, or a new one?

  # Get the current mh_proposal info.
  if ("mh_proposal" %in% names(blocks))
  {
    current_mh_proposal_info = blocks[["mh_proposal"]][[mh_proposal_number]]

    if ("mh_proposal" %in% names(current_mh_proposal_info))
    {
      # mh_proposal is complete
      print_mh_proposal_info(mh_proposal_number,blocks,line_counter-1)
      mh_proposal_number = mh_proposal_number + 1
    }
    else if ( ("evaluate_log_mh_proposal" %in% names(current_mh_proposal_info)) || ("simulate_mh_proposal" %in% names(current_mh_proposal_info)) || ("mh_transform" %in% names(current_mh_proposal_info)) || ("mh_inverse_transform" %in% names(current_mh_proposal_info)) || ("mh_transform_jacobian_matrix" %in% names(current_mh_proposal_info)) )
    {
      if ( ("mh_transform" %in% names(current_mh_proposal_info)) || ("mh_inverse_transform" %in% names(current_mh_proposal_info)) || ("mh_transform_jacobian_matrix" %in% names(current_mh_proposal_info)) )
      {
        # mh_proposal is complete
        if ( ("evaluate_log_mh_proposal" %in% names(current_mh_proposal_info)) && ("simulate_mh_proposal" %in% names(current_mh_proposal_info)) && ("mh_transform" %in% names(current_mh_proposal_info)) && ("mh_inverse_transform" %in% names(current_mh_proposal_info)) && ("mh_transform_jacobian_matrix" %in% names(current_mh_proposal_info)) )
        {
          print_mh_proposal_info(mh_proposal_number,blocks,line_counter-1)
          mh_proposal_number = mh_proposal_number + 1
        }
      }
      else
      {
        # mh_proposal is complete
        if ( ("evaluate_log_mh_proposal" %in% names(current_mh_proposal_info)) && ("simulate_mh_proposal" %in% names(current_mh_proposal_info)) )
        {
          print_mh_proposal_info(mh_proposal_number,blocks,line_counter-1)
          mh_proposal_number = mh_proposal_number + 1
        }
      }

      if (block_name=="mh_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', mh_proposal specified, but previous importance proposal block is incomplete (did you specify both evaluate_log_mh_proposal and simulate_mh_proposal, and if proposing on a transformed space, the transform, inverse transform and jacobian_matrix?)',sep=""))
      }

    }
  }
  else
  {
    mh_proposal_number = mh_proposal_number + 1
  }

  return(mh_proposal_number)
}

independent_mh_proposal_processing = function(independent_mh_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current independent_mh_proposal, or a new one?

  # Get the current independent_mh_proposal info.
  if ("independent_mh_proposal" %in% names(blocks))
  {
    current_independent_mh_proposal_info = blocks[["independent_mh_proposal"]][[independent_mh_proposal_number]]

    if ("independent_mh_proposal" %in% names(current_independent_mh_proposal_info))
    {
      # independent_mh_proposal is complete
      print_independent_mh_proposal_info(independent_mh_proposal_number,blocks,line_counter-1)
      independent_mh_proposal_number = independent_mh_proposal_number + 1
    }
    else if ( ("evaluate_log_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) || ("simulate_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) || ("independent_mh_transform" %in% names(current_mh_independent_proposal_info)) || ("independent_mh_inverse_transform" %in% names(current_mh_independent_proposal_info)) || ("independent_mh_transform_jacobian_matrix" %in% names(current_mh_independent_proposal_info)) )
    {
      if ( ("independent_mh_transform" %in% names(current_independent_mh_proposal_info)) || ("independent_mh_inverse_transform" %in% names(current_independent_mh_proposal_info)) || ("independent_mh_transform_jacobian_matrix" %in% names(current_independent_mh_proposal_info)) )
      {
        # independent_mh_proposal is complete
        if ( ("evaluate_log_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) && ("simulate_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) && ("independent_mh_transform" %in% names(current_independent_mh_proposal_info)) && ("independent_mh_inverse_transform" %in% names(current_independent_mh_proposal_info)) && ("independent_mh_transform_jacobian_matrix" %in% names(current_independent_mh_proposal_info)) )
        {
          print_independent_mh_proposal_info(independent_mh_proposal_number,blocks,line_counter-1)
          independent_mh_proposal_number = independent_mh_proposal_number + 1
        }
      }
      else
      {
        # independent_mh_proposal is complete
        if ( ("evaluate_log_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) && ("simulate_independent_mh_proposal" %in% names(current_independent_mh_proposal_info)) )
        {
          print_independent_mh_proposal_info(independent_mh_proposal_number,blocks,line_counter-1)
          independent_mh_proposal_number = independent_mh_proposal_number + 1
        }
      }

      if (block_name=="independent_mh_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', independent_mh_proposal specified, but previous importance proposal block is incomplete (did you specify both evaluate_log_independent_mh_proposal and simulate_independent_mh_proposaland if proposing on a transformed space, the transform, inverse transform and jacobian_matrix?)',sep=""))
      }

    }
  }
  else
  {
    independent_mh_proposal_number = independent_mh_proposal_number + 1
  }

  return(independent_mh_proposal_number)
}

m_proposal_processing = function(m_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current m_proposal, or a new one?

  # Get the current m_proposal info.
  if ("m_proposal" %in% names(blocks))
  {
    current_m_proposal_info = blocks[["m_proposal"]][[m_proposal_number]]

    if ("m_proposal" %in% names(current_m_proposal_info))
    {
      # m_proposal is complete
      print_m_proposal_info(m_proposal_number,blocks,line_counter-1)
      m_proposal_number = m_proposal_number + 1
    }
    else if ( ("simulate_m_proposal" %in% names(current_m_proposal_info)) || ("m_transform" %in% names(current_m_proposal_info)) || ("m_inverse_transform" %in% names(current_m_proposal_info)) || ("m_transform_jacobian_matrix" %in% names(current_m_proposal_info)) )
    {
      if ( ("m_transform" %in% names(current_m_proposal_info)) || ("m_inverse_transform" %in% names(current_m_proposal_info)) || ("m_transform_jacobian_matrix" %in% names(current_m_proposal_info)) )
      {
        # m_proposal is complete
        if ( ("simulate_m_proposal" %in% names(current_m_proposal_info)) && ("m_transform" %in% names(current_m_proposal_info)) && ("m_inverse_transform" %in% names(current_m_proposal_info)) && ("m_transform_jacobian_matrix" %in% names(current_m_proposal_info)) )
        {
          print_m_proposal_info(m_proposal_number,blocks,line_counter-1)
          m_proposal_number = m_proposal_number + 1
        }
      }
      else
      {
        # m_proposal is complete
        if ( ("simulate_m_proposal" %in% names(current_m_proposal_info)) )
        {
          print_m_proposal_info(m_proposal_number,blocks,line_counter-1)
          m_proposal_number = m_proposal_number + 1
        }
      }

      if (block_name=="m_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', m_proposal specified, but previous importance proposal block is incomplete (did you specify simulate_m_proposal,and if proposing on a transformed space, the transform, inverse transform and jacobian_matrix?)',sep=""))
      }

    }
  }
  else
  {
    m_proposal_number = m_proposal_number + 1
  }

  return(m_proposal_number)
}

transition_proposal_processing = function(transition_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current transition_proposal, or a new one?

  # Get the current transition_proposal info.
  if ("transition_proposal" %in% names(blocks))
  {
    current_transition_proposal_info = blocks[["transition_proposal"]][[transition_proposal_number]]

    if ("transition_proposal" %in% names(current_transition_proposal_info))
    {
      # transition_proposal is complete
      print_transition_proposal_info(transition_proposal_number,blocks,line_counter-1)
      transition_proposal_number = transition_proposal_number + 1
    }
    else if ( ("evaluate_log_transition_proposal" %in% names(current_transition_proposal_info)) || ("simulate_transition_proposal" %in% names(current_transition_proposal_info)) )
    {
      # transition_proposal is complete
      if ( ("evaluate_log_transition_proposal" %in% names(current_transition_proposal_info)) && ("simulate_transition_proposal" %in% names(current_transition_proposal_info)) )
      {
        print_transition_proposal_info(transition_proposal_number,blocks,line_counter-1)
        transition_proposal_number = transition_proposal_number + 1
      }

      if (block_name=="transition_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', transition_proposal specified, but previous importance proposal block is incomplete (did you specify both evaluate_log_transition_proposal and simulate_transition_proposal?)',sep=""))
      }

    }
  }
  else
  {
    transition_proposal_number = transition_proposal_number + 1
  }

  return(transition_proposal_number)
}

enk_transform_processing = function(enk_transform_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current enk_transform, or a new one?

  # Get the current enk_transform info.
  if ("enk_transform" %in% names(blocks))
  {
    current_enk_transform_info = blocks[["enk_transform"]][[enk_transform_number]]

    if ("enk_transform" %in% names(current_enk_transform_info))
    {
      # enk_transform is complete
      print_enk_transform_info(enk_transform_number,blocks,line_counter-1)
      enk_transform_number = enk_transform_number + 1
    }
    else if ( ("enk_transform" %in% names(current_enk_transform_info)) || ("enk_inverse_transform" %in% names(current_enk_transform_info)) )
    {
      # enk_transform is complete
      if ( ("enk_transform" %in% names(current_enk_transform_info)) && ("enk_inverse_transform" %in% names(current_enk_transform_info)) )
      {
        print_enk_transform_info(enk_transform_number,blocks,line_counter-1)
        enk_transform_number = enk_transform_number + 1
      }

      if (block_name=="enk_transform")
      {
        stop(paste("Invalid file: line ",line_counter,', enk_transform specified, but previous importance proposal block is incomplete (did you specify both enk_transform and enk_inverse_transform?)',sep=""))
      }

    }
  }
  else
  {
    enk_transform_number = enk_transform_number + 1
  }

  return(enk_transform_number)
}

reinforce_gradient_processing = function(reinforce_gradient_number,blocks,block_name,line_counter)
{
  # if (reinforce_gradient_number!=0)
  # {
  #   stop(paste("Invalid file: line ",line_counter,', reinforce_gradient already specified.',sep=""))
  # }
  if (reinforce_gradient_number>0)
  {
    print_reinforce_gradient_info(reinforce_gradient_number,blocks,line_counter-1)
  }
  reinforce_gradient_number = reinforce_gradient_number + 1
  return(reinforce_gradient_number)
}

method_processing = function(method_number,blocks,block_name,line_counter)
{
  # if (method_number!=0)
  # {
  #   stop(paste("Invalid file: line ",line_counter,', method already specified.',sep=""))
  # }
  if (method_number>0)
  {
    print_method_info(method_number,blocks,line_counter-1)
  }
  method_number = method_number + 1
  return(method_number)
}

print_factor_info = function(factor_index,blocks,line_counter)
{
  factor_info_string = ""
  last_factor_names = names(blocks[["factor"]][[factor_index]])
  for (j in 1:length(last_factor_names))
  {
    factor_info_string = paste(factor_info_string,last_factor_names[j],sep=", ")
  }
  if (nchar(factor_info_string)>2)
  {
    factor_info_string = substr(factor_info_string,3,nchar(factor_info_string))
  }
  print(paste('Factor ends on line ',line_counter,'. Contains ',factor_info_string,'.',sep = ""))
}

print_transition_model_info = function(transition_model_index,blocks,line_counter)
{
  transition_model_info_string = ""
  last_transition_model_names = names(blocks[["transition_model"]][[transition_model_index]])
  for (j in 1:length(last_transition_model_names))
  {
    transition_model_info_string = paste(transition_model_info_string,last_transition_model_names[j],sep=", ")
  }
  if (nchar(transition_model_info_string)>2)
  {
    transition_model_info_string = substr(transition_model_info_string,3,nchar(transition_model_info_string))
  }
  print(paste('Transition model ends on line ',line_counter,'. Contains ',transition_model_info_string,'.',sep = ""))
}

print_potential_function_info = function(potential_function_index,blocks,line_counter)
{
  potential_function_info_string = ""
  last_potential_function_names = names(blocks[["potential_function"]][[potential_function_index]])
  for (j in 1:length(last_potential_function_names))
  {
    potential_function_info_string = paste(potential_function_info_string,last_potential_function_names[j],sep=", ")
  }
  if (nchar(potential_function_info_string)>2)
  {
    potential_function_info_string = substr(potential_function_info_string,3,nchar(potential_function_info_string))
  }
  print(paste('Potential function ends on line ',line_counter,'. Contains ',potential_function_info_string,'.',sep = ""))
}

print_importance_proposal_info = function(importance_proposal_index,blocks,line_counter)
{
  importance_proposal_info_string = ""
  last_importance_proposal_names = names(blocks[["importance_proposal"]][[importance_proposal_index]])
  for (j in 1:length(last_importance_proposal_names))
  {
    importance_proposal_info_string = paste(importance_proposal_info_string,last_importance_proposal_names[j],sep=", ")
  }
  if (nchar(importance_proposal_info_string)>2)
  {
    importance_proposal_info_string = substr(importance_proposal_info_string,3,nchar(importance_proposal_info_string))
  }
  print(paste('Importance_proposal ends on line ',line_counter,'. Contains ',importance_proposal_info_string,'.',sep = ""))
}

print_mh_proposal_info = function(mh_proposal_index,blocks,line_counter)
{
  mh_proposal_info_string = ""
  last_mh_proposal_names = names(blocks[["mh_proposal"]][[mh_proposal_index]])
  for (j in 1:length(last_mh_proposal_names))
  {
    mh_proposal_info_string = paste(mh_proposal_info_string,last_mh_proposal_names[j],sep=", ")
  }
  if (nchar(mh_proposal_info_string)>2)
  {
    mh_proposal_info_string = substr(mh_proposal_info_string,3,nchar(mh_proposal_info_string))
  }
  print(paste('mh_proposal ends on line ',line_counter,'. Contains ',mh_proposal_info_string,'.',sep = ""))
}

print_independent_mh_proposal_info = function(independent_mh_proposal_index,blocks,line_counter)
{
  independent_mh_proposal_info_string = ""
  last_independent_mh_proposal_names = names(blocks[["independent_mh_proposal"]][[independent_mh_proposal_index]])
  for (j in 1:length(last_independent_mh_proposal_names))
  {
    independent_mh_proposal_info_string = paste(independent_mh_proposal_info_string,last_independent_mh_proposal_names[j],sep=", ")
  }
  if (nchar(independent_mh_proposal_info_string)>2)
  {
    independent_mh_proposal_info_string = substr(independent_mh_proposal_info_string,3,nchar(independent_mh_proposal_info_string))
  }
  print(paste('independent_mh_proposal ends on line ',line_counter,'. Contains ',independent_mh_proposal_info_string,'.',sep = ""))
}

print_m_proposal_info = function(m_proposal_index,blocks,line_counter)
{
  m_proposal_info_string = ""
  last_m_proposal_names = names(blocks[["m_proposal"]][[m_proposal_index]])
  for (j in 1:length(last_m_proposal_names))
  {
    m_proposal_info_string = paste(m_proposal_info_string,last_m_proposal_names[j],sep=", ")
  }
  if (nchar(m_proposal_info_string)>2)
  {
    m_proposal_info_string = substr(m_proposal_info_string,3,nchar(m_proposal_info_string))
  }
  print(paste('m_proposal ends on line ',line_counter,'. Contains ',m_proposal_info_string,'.',sep = ""))
}

print_transition_proposal_info = function(transition_proposal_index,blocks,line_counter)
{
  transition_proposal_info_string = ""
  last_transition_proposal_names = names(blocks[["transition_proposal"]][[transition_proposal_index]])
  for (j in 1:length(last_transition_proposal_names))
  {
    transition_proposal_info_string = paste(transition_proposal_info_string,last_transition_proposal_names[j],sep=", ")
  }
  if (nchar(transition_proposal_info_string)>2)
  {
    transition_proposal_info_string = substr(transition_proposal_info_string,3,nchar(transition_proposal_info_string))
  }
  print(paste('transition_proposal ends on line ',line_counter,'. Contains ',transition_proposal_info_string,'.',sep = ""))
}

print_enk_transform_info = function(enk_transform_index,blocks,line_counter)
{
  enk_transform_info_string = ""
  last_enk_transform_names = names(blocks[["enk_transform"]][[enk_transform_index]])
  for (j in 1:length(last_enk_transform_names))
  {
    enk_transform_info_string = paste(enk_transform_info_string,last_enk_transform_names[j],sep=", ")
  }
  if (nchar(enk_transform_info_string)>2)
  {
    enk_transform_info_string = substr(enk_transform_info_string,3,nchar(enk_transform_info_string))
  }
  print(paste('enk_transform ends on line ',line_counter,'. Contains ',enk_transform_info_string,'.',sep = ""))
}

print_reinforce_gradient_info = function(reinforce_gradient_index,blocks,line_counter)
{
  reinforce_gradient_info_string = ""
  last_reinforce_gradient_names = names(blocks[["reinforce_gradient"]][[reinforce_gradient_index]])
  for (j in 1:length(last_reinforce_gradient_names))
  {
    reinforce_gradient_info_string = paste(reinforce_gradient_info_string,last_reinforce_gradient_names[j],sep=", ")
  }
  if (nchar(reinforce_gradient_info_string)>2)
  {
    reinforce_gradient_info_string = substr(reinforce_gradient_info_string,3,nchar(reinforce_gradient_info_string))
  }
  print(paste('reinforce_gradient ends on line ',line_counter,'. Contains ',reinforce_gradient_info_string,'.',sep = ""))
}

print_method_info = function(method_index,blocks,line_counter)
{
  method_info_string = ""
  last_method_names = names(blocks[["method"]][[method_index]])
  for (j in 1:length(last_method_names))
  {
    method_info_string = paste(method_info_string,last_method_names[j],sep=", ")
  }
  if (nchar(method_info_string)>2)
  {
    method_info_string = substr(method_info_string,3,nchar(method_info_string))
  }
  print(paste('Method ends on line ',line_counter,'. Contains ',method_info_string,'.',sep = ""))
}

print_data_info = function(data_index,blocks,line_counter)
{
  data_info_string = ""
  last_data_names = names(blocks[["data"]][[data_index]])
  for (j in 1:length(last_data_names))
  {
    data_info_string = paste(data_info_string,last_data_names[j],sep=", ")
  }
  if (nchar(data_info_string)>2)
  {
    data_info_string = substr(data_info_string,3,nchar(data_info_string))
  }
  print(paste('data ends on line ',line_counter,'. Contains ',data_info_string,'.',sep = ""))
}

determine_block_type = function(split_block_name,blocks,line_counter,block_type,block_name,factor_number,transition_model_number,potential_function_number,importance_proposal_number,mh_proposal_number,independent_mh_proposal_number,m_proposal_number,enk_transform_number,transition_proposal_number,data_number,method_number)
{
  if (length(split_block_name)==1)
  {
    is_custom = TRUE
    block_name = split_block_name
    block_function = ""
  }
  else if (length(split_block_name)==2)
  {
    is_custom = FALSE
    block_name = split_block_name[1]
    block_function = split_block_name[2]
  }
  else
  {
    stop(paste("Invalid file: line ",line_counter,", invalid function definition (can only have up two parts, separated by commas).",sep=""))
  }

  prior_function_types = c("prior","evaluate_log_prior","simulate_prior","evaluate_gradient_log_prior","evaluate_second_gradient_log_prior")
  custom_likelihood_function_types = c("evaluate_log_likelihood","evaluate_gradient_log_likelihood","evaluate_second_gradient_log_likelihood")
  #file_likelihood_types = c("importance_sample","smc_mcmc_move")
  sbi_likelihood_function_types = c("simulate_data_model","sbi_likelihood","summary_statistics")
  linear_gaussian_data_model_types = c("linear_gaussian_data_model","linear_gaussian_data_matrix","linear_gaussian_data_covariance")
  nonlinear_gaussian_data_model_types = c("nonlinear_gaussian_data_model","nonlinear_gaussian_data_function","nonlinear_gaussian_data_covariance")
  other_likelihood_function_types = c("likelihood")
  factor_function_types = c(prior_function_types,custom_likelihood_function_types,sbi_likelihood_function_types,linear_gaussian_data_model_types,nonlinear_gaussian_data_model_types,other_likelihood_function_types)
  ilike_transition_model_types = c("transition_model")
  linear_gaussian_transition_model_types = c("linear_gaussian_transition_model","linear_gaussian_transition_matrix","linear_gaussian_transition_covariance")
  nonlinear_gaussian_transition_model_types = c("nonlinear_gaussian_transition_model","nonlinear_gaussian_transition_function","nonlinear_gaussian_transition_covariance")
  custom_transition_model_types = c("simulate_transition_model","evaluate_log_transition_model")
  transition_model_types = c(ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_types)
  custom_potential_function_types = c("evaluate_log_potential_function")
  ilike_potential_function_types = c("potential_function")
  potential_function_types = c(custom_potential_function_types,ilike_potential_function_types)
  data_function_types = c("data")
  importance_proposal_types = c("simulate_importance_proposal","evaluate_log_importance_proposal")
  mh_proposal_types = c("simulate_mh_proposal","evaluate_log_mh_proposal","mh_proposal","mh_transform","mh_inverse_transform","mh_transform_jacobian_matrix")
  independent_mh_proposal_types = c("simulate_independent_mh_proposal","evaluate_log_independent_mh_proposal","independent_mh_proposal","independent_mh_transform","independent_mh_inverse_transform","independent_mh_transform_jacobian_matrix")
  m_proposal_types = c("simulate_m_proposal","m_proposal","m_transform","m_inverse_transform","m_transform_jacobian_matrix")
  transition_proposal_types = c("simulate_transition_proposal","evaluate_log_transition_proposal")
  enk_transform_types = c("enk_transform","enk_inverse_transform")
  #reinforce_gradient = c("reinforce_gradient")
  method_function_types = c("mcmc_weights","mcmc_termination","adaptive_resampling","adaptive_target","smc_termination","smc_sequence","reinforce_gradient")

  # distinguish between factor, data, etc
  if (block_name %in% factor_function_types)
  {
    block_type = "factor"
    factor_number = factor_processing(factor_number,blocks,block_name,prior_function_types,custom_likelihood_function_types,sbi_likelihood_function_types,linear_gaussian_data_model_types,nonlinear_gaussian_data_model_types,other_likelihood_function_types,line_counter)
    number_to_pass_to_extract_block = factor_number
  }
  else if (block_name %in% transition_model_types)
  {
    block_type = "transition_model"
    transition_model_number = transition_model_processing(transition_model_number,blocks,block_name,ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_types,line_counter)
    number_to_pass_to_extract_block = transition_model_number
  }
  else if (block_name %in% potential_function_types)
  {
    block_type = "potential_function"
    potential_function_number = potential_function_processing(potential_function_number,blocks,block_name,custom_potential_function_types,ilike_potential_function_types,line_counter)
    number_to_pass_to_extract_block = potential_function_number
  }
  else if (block_name %in% data_function_types)
  {
    block_type = "data"
    data_number = data_processing(data_number,line_counter)
    number_to_pass_to_extract_block = data_number
  }
  else if (block_name %in% importance_proposal_types)
  {
    block_type = "importance_proposal"
    importance_proposal_number = importance_proposal_processing(importance_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = importance_proposal_number
  }
  else if (block_name %in% mh_proposal_types)
  {
    block_type = "mh_proposal"
    mh_proposal_number = mh_proposal_processing(mh_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = mh_proposal_number
  }
  else if (block_name %in% independent_mh_proposal_types)
  {
    block_type = "independent_mh_proposal"
    independent_mh_proposal_number = independent_mh_proposal_processing(independent_mh_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = independent_mh_proposal_number
  }
  else if (block_name %in% m_proposal_types)
  {
    block_type = "m_proposal"
    m_proposal_number = m_proposal_processing(m_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = m_proposal_number
  }
  else if (block_name %in% enk_transform_types)
  {
    block_type = "enk_transform"
    enk_transform_number = enk_transform_processing(enk_transform_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = enk_transform_number
  }
  else if (block_name %in% transition_proposal_types)
  {
    block_type = "transition_proposal"
    transition_proposal_number = transition_proposal_processing(transition_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = transition_proposal_number
  }
  #else if (block_name %in% reinforce_gradient_types)
  #{
  #  block_type = "reinforce_gradient"
  #  reinforce_gradient_number = reinforce_gradient_processing(reinforce_gradient_number,blocks,block_name,line_counter)
  #  number_to_pass_to_extract_block = reinforce_gradient_number
  #}
  else if (block_name %in% method_function_types)
  {
    block_type = "method"
    method_number = method_processing(method_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = method_number
  }
  else
  {
    stop(paste("Invalid file: line ",line_counter,', function type "',block_name,'" unknown.',sep=""))
  }

  return(list(block_type,
              number_to_pass_to_extract_block,
              block_name,
              is_custom,
              block_function,
              factor_number,
              transition_model_number,
              potential_function_number,
              importance_proposal_number,
              mh_proposal_number,
              independent_mh_proposal_number,
              m_proposal_number,
              enk_transform_number,
              transition_proposal_number,
              data_number,
              method_number))
}

extract_block <- function(blocks,block_type,block_name,factor_number,line_counter,block_code,block_function,is_custom,parameter_list,external_packages,julia_bin_dir,julia_required_libraries)
{
  # Get information about the order in which MCMC moves are included.
  if ( (block_type=="mh_proposal") || (block_type=="independent_mh_proposal") || (block_type=="m_proposal") )
  {
    # Are we adding to an existing block? If we are not, then we need to add in the information about the order of the new MCMC move.
    if (factor_number>length(blocks[[block_type]]))
    {
      if (block_type=="mh_proposal")
      {
        blocks[["order_of_mcmc"]] = c(blocks[["order_of_mcmc"]],1)
      }

      if (block_type=="independent_mh_proposal")
      {
        blocks[["order_of_mcmc"]] = c(blocks[["order_of_mcmc"]],2)
      }

      if (block_type=="m_proposal")
      {
        blocks[["order_of_mcmc"]] = c(blocks[["order_of_mcmc"]],3)
      }
    }
  }

  # ignore block if block number is not positive
  if (factor_number>0)
  {
    if (is_custom==TRUE)
    {
      my_list = list(RcppXPtrUtils::cppXPtr(block_code,plugins=c("cpp11"),depends = paste(c("ilike","RcppArmadillo","BH","dqrng","sitmo"),external_packages) ))
      names(my_list) = c(block_name)

      if (length((blocks[[block_type]]))==0)
      {
        blocks[[block_type]][[factor_number]] = my_list
      }
      else if (factor_number!=length((blocks[[block_type]])))
      {
        blocks[[block_type]][[factor_number]] = my_list
      }
      else
      {
        blocks[[block_type]][[factor_number]] = append(blocks[[block_type]][[factor_number]],my_list)
      }

    }
    else
    {
      # Check if block_function is of the form: function(some,stuff)

      block_function <- gsub("\\s+", "", block_function)

      # Split input into function name and arguments
      parts <- strsplit(block_function, "\\(")[[1]]
      if (length(parts)>1)
      {
        function_name <- parts[1]
        if ( (substr(block_function,nchar(function_name)+1,nchar(function_name)+1)!="(") || (substr(block_function,nchar(block_function),nchar(block_function))!=")"))
        {
          stop(paste("Block ",block_name,", line number ",line_counter,": input does not specify a function.",sep=""))
        }
        arg_string <- substr(block_function,nchar(function_name)+2,nchar(block_function)-1)

        ilike_split_name = strsplit(function_name,"::")[[1]]

        if ( (length(ilike_split_name)>1) && (ilike_split_name[1]=="ilike") )
        {
          is_like_function = TRUE
          ilike_type = paste(ilike_split_name[2:length(ilike_split_name)],collapse = "::")
        }
        else
        {
          is_like_function = FALSE
        }
      }
      else
      {
        stop(paste("Block ",block_name,", line number ",line_counter,": function needs to have a name.",sep=""))
      }

      if (is_like_function==TRUE)
      {
        if ( (block_name=="prior") || (block_name=="importance_proposal") || (block_name=="independent_mh_proposal") || (block_name=="mh_proposal") || (block_name=="m_proposal") || (block_name=="sbi_likelihood") || (block_name=="smc_sequence") || (block_name=="linear_gaussian_transition_model") || (block_name=="nonlinear_gaussian_transition_model") || (block_name=="linear_gaussian_data_model") || (block_name=="nonlinear_gaussian_data_model") )
        {
          split_arg_string = split_string_at_comma_ignoring_parentheses(arg_string)

          variables = split_arg_string[1]

          parameters = list()
          if (length(split_arg_string)>1)
          {
            for (k in 2:(length(split_arg_string)))
            {
              parameters[[k-1]] = split_arg_string[k]
            }
          }

          my_list = list(list(type=ilike_type,
                              variables=variables,
                              parameters=parameters))
          names(my_list) = c(block_name)

          if (length((blocks[[block_type]]))==0)
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else if (factor_number!=length((blocks[[block_type]])))
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else
          {
            blocks[[block_type]][[factor_number]] = append(blocks[[block_type]][[factor_number]],my_list)
          }

        }
        else if ( (block_name=="mcmc_termination") || (block_name=="adaptive_resampling") || (block_name=="adaptive_target") || (block_name=="smc_termination") || (block_name=="reinforce_gradient") )
        {
          split_arg_string = split_string_at_comma_ignoring_parentheses(arg_string)

          parameters = list()
          if (length(split_arg_string)>0)
          {
            for (k in 1:(length(split_arg_string)))
            {
              parameters[[k]] = split_arg_string[k]
            }
          }

          my_list = list(list(type=ilike_type,
                              parameters=parameters))
          names(my_list) = c(block_name)

          if (length((blocks[[block_type]]))==0)
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else if (factor_number!=length((blocks[[block_type]])))
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else
          {
            blocks[[block_type]][[factor_number]] = append(blocks[[block_type]][[factor_number]],my_list)
          }

        }
        else if (block_name=="likelihood")
        {
          split_arg_string = split_string_at_comma_ignoring_parentheses(arg_string)

          filename = split_arg_string[1]

          parameters = list()
          if (length(split_arg_string)>1)
          {
            for (k in 2:(length(split_arg_string)))
            {
              parameters[[k-1]] = split_arg_string[k]
            }
          }

          my_list = list(list(type=ilike_type,
                              model=ilike::parse_ilike_model(filename,parameter_list,external_packages,julia_bin_dir,julia_required_libraries),
                              parameters=parameters))
          names(my_list) = c(block_name)

          if (length((blocks[[block_type]]))==0)
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else if (factor_number!=length((blocks[[block_type]])))
          {
            blocks[[block_type]][[factor_number]] = my_list
          }
          else
          {
            blocks[[block_type]][[factor_number]] = append(blocks[[block_type]][[factor_number]],my_list)
          }
        }
        else
        {
          stop(paste("Block ",block_name,", line number ",line_counter,": ",block_name," is invalid block type.",sep=""))
        }
      }
      else
      {

        # need to check if output is assigned from function (for simulate)
        split_at_equals = strsplit(block_function,"=")[[1]]
        if (length(split_at_equals)==2)
        {
          output_variable = split_at_equals[1]
          function_to_use = split_at_equals[2]
          function_name = strsplit(function_name,"=")[[1]][2]
        }
        else
        {
          output_variable = ""
          function_to_use = block_function
        }

        if (length(split_at_equals)>2)
        {
          stop(paste("Block ",block_name,", line number ",line_counter,": can have a maximum of one equals sign in it.",sep=""))
        }

        function_info = ilike_parse(function_to_use,
                                    parameter_list)

        R_function_arguments = function_info[[2]]

        cpp_function_arguments_string = ""
        which_proposed_parameters = c()
        which_parameters = c()
        which_proposal_parameters = c()
        which_data = c()
        if (length(R_function_arguments)>0)
        {
          which_proposed_parameters = matrix(0,length(R_function_arguments))
          which_parameters = matrix(0,length(R_function_arguments))
          which_proposal_parameters = matrix(0,length(R_function_arguments))
          which_data = matrix(0,length(R_function_arguments))

          for (i in 1:length(R_function_arguments))
          {
            split_at_dot = strsplit(R_function_arguments[i],"\\.")[[1]]
            if (length(split_at_dot)==0)
            {
              stop(paste("Block ",block_name,", line number ",line_counter,": argument is of size zero.",sep=""))
            }
            else if (length(split_at_dot)==1)
            {
              cpp_function_arguments_string = paste(cpp_function_arguments_string,split_at_dot[1],sep=",")
            }
            else if (length(split_at_dot)>=2)
            {
              cpp_function_arguments_string = paste(cpp_function_arguments_string,',',split_at_dot[1],'["',paste(split_at_dot[2:length(split_at_dot)],collapse=""),'"]',sep="")

              if (split_at_dot[1]=="proposed_parameters")
              {
                which_proposed_parameters[i] = 1
              }

              if (split_at_dot[1]=="parameters")
              {
                which_parameters[i] = 1
              }

              if (split_at_dot[1]=="proposal_parameters")
              {
                which_proposal_parameters[i] = 1
              }

              if (split_at_dot[1]=="data")
              {
                which_data[i] = 1
              }
            }
          }

          cpp_function_arguments_string = substr(cpp_function_arguments_string,2,nchar(cpp_function_arguments_string))
        }

        # give C++ wrapper for this
        # write to cpp file
        # source
        # delete file

        cpp_function_name = paste(block_name,factor_number,sep="")
        R_function_name = paste(cpp_function_name,'_R',sep="")

        temp_filename = paste(cpp_function_name,".cpp",sep="")
        xptr_name = paste(cpp_function_name,"getXPtr",sep="_")

        proposal_type = 0

        if (block_name=="data")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a data function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in a data function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in a data function.",sep=""))
          }

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a prior.",sep=""))
          }

          return_type = "Data"
          arguments = c()
          #args_for_typedef = ""
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Data output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="evaluate_log_prior")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a prior.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="evaluate_log_likelihood")
        {
          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="simulate_prior")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a prior.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="evaluate_log_importance_proposal")
        {

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0)  )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            arguments[3] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 1
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="simulate_importance_proposal")
        {
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()

          if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Parameters &proposal_parameters"
            arguments[3] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 2
          }
          else
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            #args_for_typedef = "const Parameters&"

            proposal_type = 1
          }

          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="evaluate_log_mh_proposal")
        {
          if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &parameters"
            arguments[3] = "const Parameters &proposal_parameters"
            arguments[4] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &parameters"
            arguments[3] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &parameters"
            arguments[3] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 1
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="simulate_mh_proposal")
        {
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()

          if (sum(which_parameters)>0)
          {
            if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              arguments[4] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"

              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              #args_for_typedef = "const Parameters&"

              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, should this have been specified as simulate_independent_mh_proposal?",sep=""))
          }

          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="evaluate_log_independent_mh_proposal")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            arguments[3] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"

            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            #args_for_typedef = "const Parameters&"

            proposal_type = 1
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="simulate_independent_mh_proposal")
        {
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()

          if (sum(which_parameters)==0)
          {
            if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &proposal_parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"

              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              #args_for_typedef = "const Parameters&"

              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": parameters used, should this have been specified as simulate_mh_proposal?",sep=""))
          }

          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="simulate_m_proposal")
        {
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()

          if (sum(which_parameters)>0)
          {
            if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              arguments[4] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"

              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"

              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              #args_for_typedef = "const Parameters&"

              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, invalid proposal type?",sep=""))
          }

          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")

        }
        else if (block_name=="mcmc_weights")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a mcmc_weights function.",sep=""))
          }

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a mcmc_weights function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "NumericVector"
          arguments = c()

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); NumericVector output = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="simulate_data_model")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a data_model.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "Data"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          arguments[2] = "const Parameters &parameters"
          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Data output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="summary_statistics")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Data"
          arguments = c()
          arguments[1] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Data output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="nonlinear_gaussian_data_function")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Data"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Data output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="simulate_transition_model")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          arguments[2] = "const Parameters &parameters"
          arguments[3] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="nonlinear_gaussian_transition_function")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="simulate_transition_proposal")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          arguments[2] = "const Parameters &parameters"
          arguments[3] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="enk_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="enk_inverse_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="linear_gaussian_data_matrix")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="linear_gaussian_data_covariance")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="nonlinear_gaussian_data_covariance")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="linear_gaussian_transition_matrix")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="linear_gaussian_transition_covariance")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="nonlinear_gaussian_transition_covariance")
        {

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="evaluate_log_transition_model")
        {

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return NumericVector(f(',cpp_function_arguments_string,'))[0]; return output;',sep="")
        }
        else if (block_name=="evaluate_log_potential_function")
        {

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return NumericVector(f(',cpp_function_arguments_string,'))[0]; return output;',sep="")
        }
        else if (block_name=="evaluate_log_transition_proposal")
        {

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          arguments[2] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return NumericVector(f(',cpp_function_arguments_string,'))[0]; return output;',sep="")
        }
        else if (block_name=="m_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="m_inverse_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="m_transform_jacobian_matrix")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="mh_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="mh_inverse_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="mh_transform_jacobian_matrix")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else if (block_name=="independent_mh_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="independent_mh_inverse_transform")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="independent_mh_transform_jacobian_ma")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_proposed_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposed parameters not possible in this function.",sep=""))
          }

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "arma::mat"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return f(',cpp_function_arguments_string,');',sep="")
        }
        else
        {
          stop(paste("Block ",block_name,", line number ",line_counter,": block is of unknown type.",sep=""))
        }

        args_for_typedef = paste(arguments,collapse=",")

        code = paste(return_type,' ',cpp_function_name,'(',sep="")
        if (length(arguments)>0)
        {
          for (i in 1:length(arguments))
          {
            code = paste(code,arguments[i],sep="")
            if (i!=length(arguments))
            {
              code = paste(code,',',sep="")
            }
          }
        }
        code = paste(code,') {',function_body,' }',sep="")

        fileConn<-file(temp_filename)

        writeLines(c(
          '#include <RcppArmadillo.h>',
          '// [[Rcpp::depends(RcppArmadillo)]]',
          '// [[Rcpp::depends(ilike)]]',
          '// [[Rcpp::depends(BH)]]',
          '// [[Rcpp::depends(dqrng)]]',
          '// [[Rcpp::depends(sitmo)]]',
          '/*** R',
          paste(R_function_name,function_info[[1]],sep=""),
          '*/',
          'using namespace Rcpp;',
          '#include <ilike.h>',
          paste("SEXP ",xptr_name,"();",sep=""),
          code,
          "// [[Rcpp::export]]",
          paste("SEXP ",xptr_name,"() {",sep=""),
          paste("  typedef", return_type, "(*funcPtr)(", args_for_typedef, ");"),
          paste("  return XPtr<funcPtr>(new funcPtr(&", cpp_function_name, "));"),
          "}"), fileConn)

        # writeLines(c(
        #   '/*** R',
        #   paste(R_function_name,function_info[[1]],sep=""),
        #   '*/'), fileConn)

        close(fileConn)

        Rcpp::sourceCpp(temp_filename)
        file.remove(temp_filename)

        my_list = list(get(xptr_name)())
        #my_list = list(RcppXPtrUtils::cppXPtr(code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo")))

        names(my_list) = c(block_name)

        if (proposal_type!=0)
          my_list[["type"]] = proposal_type

        if ("type" %in% names(blocks[[block_type]][[factor_number]]))
        {
          if (blocks[[block_type]][[factor_number]][["type"]]!=proposal_type)
            stop(paste("Block ",block_name,", line number ",line_counter,": different parts of the proposal are of different types (maybe one relies on additional parameters and/or data, and the other does not?)",sep=""))
        }

        if (length((blocks[[block_type]]))==0)
        {
          blocks[[block_type]][[factor_number]] = my_list
        }
        else if (factor_number!=length((blocks[[block_type]])))
        {
          blocks[[block_type]][[factor_number]] = my_list
        }
        else
        {
          blocks[[block_type]][[factor_number]] = append(blocks[[block_type]] [[factor_number]],my_list)
        }
      }

    }
  }

  return(blocks)
}


check_types = function(blocks)
{
  if ("data" %in% names(blocks))
  {
    for (i in 1:length(blocks[["data"]]))
    {
      if (inherits(blocks[["data"]][[i]][["data"]],"XPtr"))
        RcppXPtrUtils::checkXPtr(blocks[["data"]][[i]][["data"]], "Data")
    }
  }

  for (i in 1:length(blocks[["factor"]]))
  {
    current_factor = blocks[["factor"]][[i]]

    if ("evaluate_log_prior" %in% names(current_factor))
    {
      if (inherits(current_factor[["evaluate_log_prior"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_prior"]], "double", c("const Parameters&"))
      }
    }

    if ("simulate_prior" %in% names(current_factor))
    {
      if (inherits(current_factor[["simulate_prior"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["simulate_prior"]], "Parameters", c("RandomNumberGenerator&"))
      }
    }

    if ("evaluate_log_likelihood" %in% names(current_factor))
    {
      if (inherits(current_factor[["evaluate_log_likelihood"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_likelihood"]],  "double", c("const Parameters&","const Data&"))
      }
    }

    if ("simulate_data_model" %in% names(current_factor))
    {
      if (inherits(current_factor[["simulate_data_model"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["simulate_data_model"]],  "Data", c("RandomNumberGenerator&","const Parameters&"))
      }
    }

    if ("summary_statistics" %in% names(current_factor))
    {
      if (inherits(current_factor[["summary_statistics"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["summary_statistics"]],  "Data", c("const Data&"))
      }
    }

    if ("nonlinear_gaussian_data_function" %in% names(current_factor))
    {
      if (inherits(current_factor[["nonlinear_gaussian_data_function"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["nonlinear_gaussian_data_function"]],  "Data", c("const Parameters&"))
      }
    }

    if ("linear_gaussian_data_matrix" %in% names(current_factor))
    {
      if (inherits(current_factor[["linear_gaussian_data_matrix"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["linear_gaussian_data_matrix"]],  "arma::mat", c("const Parameters&"))
      }
    }

    if ("linear_gaussian_data_covariance" %in% names(current_factor))
    {
      if (inherits(current_factor[["linear_gaussian_data_covariance"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["linear_gaussian_data_covariance"]],  "arma::mat", c("const Parameters&"))
      }
    }

    if ("nonlinear_gaussian_data_covariance" %in% names(current_factor))
    {
      if (inherits(current_factor[["nonlinear_gaussian_data_covariance"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[["nonlinear_gaussian_data_covariance"]],  "arma::mat", c("const Parameters&"))
      }
    }

  }

  for (i in 1:length(blocks[["importance_proposal"]]))
  {
    current_importance_proposal = blocks[["importance_proposal"]][[i]]

    if ("evaluate_log_importance_proposal" %in% names(current_importance_proposal))
    {
      if (inherits(current_importance_proposal[["evaluate_log_importance_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid importance proposal specified.")
        }

        blocks[["importance_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
    }

    if ("simulate_importance_proposal" %in% names(current_importance_proposal))
    {
      if (inherits(current_importance_proposal[["simulate_importance_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid importance proposal specified.")
        }

        if ("type" %in% names(blocks[["importance_proposal"]][[i]]))
        {
          if (blocks[["importance_proposal"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("evaluate_log_importance_proposal and simulate_importance_proposal are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both evaluate_log_importance_proposal and simulate_importance_proposal.")
        }
      }
    }
  }

  for (i in 1:length(blocks[["mh_proposal"]]))
  {
    current_mh_proposal = blocks[["mh_proposal"]][[i]]

    if ("evaluate_log_mh_proposal" %in% names(current_mh_proposal))
    {
      if (inherits(current_mh_proposal[["evaluate_log_mh_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["evaluate_log_mh_proposal"]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["evaluate_log_mh_proposal"]],  "double", c("const Parameters&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["evaluate_log_mh_proposal"]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["evaluate_log_mh_proposal"]],  "double", c("const Parameters&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid importance proposal specified.")
        }

        blocks[["mh_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
    }

    if ("simulate_mh_proposal" %in% names(current_mh_proposal))
    {
      if (inherits(current_mh_proposal[["simulate_mh_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["simulate_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["simulate_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["simulate_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_mh_proposal[["simulate_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid mh proposal specified.")
        }

        if ("type" %in% names(blocks[["mh_proposal"]][[i]]))
        {
          if (blocks[["mh_proposal"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("evaluate_log_mh_proposal and simulate_mh_proposal are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both evaluate_log_mh_proposal and simulate_mh_proposal.")
        }
      }
    }

  }

  for (i in 1:length(blocks[["independent_mh_proposal"]]))
  {
    current_independent_mh_proposal = blocks[["independent_mh_proposal"]][[i]]

    if ("evaluate_log_independent_mh_proposal" %in% names(current_independent_mh_proposal))
    {
      if (inherits(current_independent_mh_proposal[["evaluate_log_independent_mh_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["evaluate_log_independent_mh_proposal"]],  "double", c("const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["evaluate_log_independent_mh_proposal"]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["evaluate_log_independent_mh_proposal"]],  "double", c("const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["evaluate_log_independent_mh_proposal"]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid importance proposal specified.")
        }

        blocks[["independent_mh_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
    }

    if ("simulate_independent_mh_proposal" %in% names(current_independent_mh_proposal))
    {
      if (inherits(current_independent_mh_proposal[["simulate_independent_mh_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["simulate_independent_mh_proposal"]], "Parameters", c("RandomNumberGenerator&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["simulate_independent_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["simulate_independent_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_independent_mh_proposal[["simulate_independent_mh_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid independent mh proposal specified.")
        }

        if ("type" %in% names(blocks[["independent_mh_proposal"]][[i]]))
        {
          if (blocks[["independent_mh_proposal"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("evaluate_log_independent_mh_proposal and simulate_independent_mh_proposal are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both evaluate_log_independent_mh_proposal and simulate_independent_mh_proposal.")
        }
      }
    }

  }

  for (i in 1:length(blocks[["m_proposal"]]))
  {
    current_m_proposal = blocks[["m_proposal"]][[i]]

    if ("simulate_m_proposal" %in% names(current_m_proposal))
    {
      if (inherits(current_m_proposal[["simulate_m_proposal"]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_m_proposal[["simulate_m_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_m_proposal[["simulate_m_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_m_proposal[["simulate_m_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_m_proposal[["simulate_m_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid m proposal specified.")
        }

        blocks[["m_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
    }

  }

  for (i in 1:length(blocks[["transition_model"]]))
  {
    current_transition_model = blocks[["transition_model"]][[i]]

    if ("evaluate_log_transition_model" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["evaluate_log_transition_model"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["evaluate_log_transition_model"]],  "double", c("const Parameters&","const Parameters&","const Data&"))
      }
    }

    if ("simulate_transition_model" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["simulate_transition_model"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["simulate_transition_model"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&"))
      }
    }

    if ("nonlinear_gaussian_transition_function" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["nonlinear_gaussian_transition_function"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["nonlinear_gaussian_transition_function"]],  "Data", c("const Parameters&","const Data&"))
      }
    }

    if ("linear_gaussian_transition_matrix" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["linear_gaussian_transition_matrix"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["linear_gaussian_transition_matrix"]],  "arma::mat", c("const Parameters&","const Data&"))
      }
    }

    if ("linear_gaussian_transition_covariance" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["linear_gaussian_transition_covariance"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["linear_gaussian_transition_covariance"]],  "arma::mat", c("const Parameters&","const Data&"))
      }
    }

    if ("nonlinear_gaussian_transition_covariance" %in% names(current_transition_model))
    {
      if (inherits(current_transition_model[["nonlinear_gaussian_transition_covariance"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_model[["nonlinear_gaussian_transition_covariance"]],  "arma::mat", c("const Parameters&","const Data&"))
      }
    }

  }

  for (i in 1:length(blocks[["transition_proposal"]]))
  {
    current_transition_proposal = blocks[["transition_proposal"]][[i]]

    if ("evaluate_log_transition_proposal" %in% names(current_transition_proposal))
    {
      if (inherits(current_transition_proposal[["evaluate_log_transition_proposal"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_proposal[["evaluate_log_transition_proposal"]],  "double", c("const Parameters&","const Parameters&","const Data&"))
      }
    }

    if ("simulate_transition_proposal" %in% names(current_transition_proposal))
    {
      if (inherits(current_transition_proposal[["simulate_transition_proposal"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_transition_proposal[["simulate_transition_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&"))
      }
    }

  }

  for (i in 1:length(blocks[["enk_transform"]]))
  {
    current_enk_transform = blocks[["enk_transform"]][[i]]

    if ("enk_transform" %in% names(current_enk_transform))
    {
      if (inherits(current_enk_transform[["enk_transform"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_enk_transform[["enk_transform"]],  "Parameters", c("const Parameters&"))
      }
    }

    if ("enk_inverse_transform" %in% names(current_enk_transform))
    {
      if (inherits(current_enk_transform[["enk_inverse_transform"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_enk_transform[["enk_inverse_transform"]],  "Parameters", c("const Parameters&"))
      }
    }

  }

  for (i in 1:length(blocks[["potential_function"]]))
  {
    current_potential_function = blocks[["potential_function"]][[i]]

    if ("evaluate_log_potential_function" %in% names(current_potential_function))
    {
      if (inherits(current_potential_function[["evaluate_log_potential_function"]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_potential_function[["evaluate_log_potential_function"]],  "double", c("const Parameters&","const Parameters&","const Data&"))
      }
    }

  }

  if ("method" %in% names(blocks))
  {
    for (i in 1:length(blocks[["method"]]))
    {
      if ("mcmc_weights" %in% names(blocks[["method"]][[i]]))
      {
        if (inherits(blocks[["method"]][[i]][["mcmc_weights"]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(blocks[["method"]][[i]][["mcmc_weights"]], "NumericVector")
        }
      }
    }
  }

  return(blocks)
}


#' Parse .ilike file to give ilike model.
#'
#' @param filename The name (and path) of the .ilike file containing the model.
#' @param parameter_list (optional) A list containing parameters for the model.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @return A list containing the model details.
#' @export
parse_ilike_model <- function(filename,
                              parameter_list = list(),
                              external_packages = c(),
                              julia_bin_dir="",
                              julia_required_libraries=c())
{
  basename = tools::file_path_sans_ext(filename)

  # Check if there is a file with a .py extension.
  if (file.exists(paste(basename,".py",sep="")))
  {
    reticulate::source_python(paste(basename,".py",sep=""),envir=globalenv())
  }

  # Check if there is a file with a .jl extension.
  if (file.exists(paste(basename,".jl",sep="")))
  {
    if (julia_bin_dir!="")
    {
      JuliaCall::julia_setup(julia_bin_dir)
      #Sys.setenv(JULIA_BINDIR = julia_bin_dir)
      output = my_julia_source(paste(basename,".jl",sep=""),julia_required_libraries)
      list2env(output, .GlobalEnv)
    }
    else
    {
      stop("Julia binary directory needs to be specified to use Julia functions.")
    }
    #JuliaCall::julia_source(paste(basename,".jl",sep=""))

    # Split julia file into lots of little files, and call mystring <- read_file("my_rnorm(n).jl"); my_rnorm <- JuliaCall::julia_eval(mystring) on each one
  }

  # Check if there is a file with a .R extension.
  if (file.exists(paste(basename,".R",sep="")))
  {
    source(paste(basename,".R",sep=""))
  }

  the_file = file(filename,open="r")

  in_block = FALSE
  blocks = list(order_of_mcmc=c())
  block_code = ""
  starting_block_flag = FALSE
  is_custom = FALSE
  line_counter = 0
  add_to_block_code = TRUE
  function_names = c()
  in_factor = FALSE

  factor_number = 0
  transition_model_number = 0
  potential_function_number = 0
  importance_proposal_number = 0
  mh_proposal_number = 0
  independent_mh_proposal_number = 0
  m_proposal_number = 0
  enk_transform_number = 0
  transition_proposal_number = 0
  data_number = 0
  #reinforce_gradient_number = 0
  method_number = 0
  block_type = "none"
  block_name = "none"
  #header = ""

  while ( TRUE ) {
    line = readLines(the_file, n = 1)
    line_counter = line_counter + 1
    if ( length(line) == 0 ) {
      break
    }

    if (starting_block_flag==TRUE)
    {
      starting_block_flag = FALSE
    }

    if ( (nchar(line)>=3) && (substr(line, 1, 3)=="//#") )
    {
      blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,parameter_list,external_packages,julia_bin_dir,julia_required_libraries)
      in_block = FALSE
      starting_block_flag = TRUE
    }
    else if (nchar(line)>=8)
    {
      if (substr(line, 1, 4)=="/***")
      {
        if (substr(line, nchar(line) - 4 + 1, nchar(line))=="***/")
        {
          # legitimate new section
          starting_block_flag = TRUE

          if (in_block==FALSE) # first block, or previous block wrapped up using //#
          {
            in_block = TRUE
          }
          else # end current block
          {
            blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,parameter_list,external_packages,julia_bin_dir,julia_required_libraries)
            is_custom = FALSE
            block_code = ""
          }

          unparsed_block_name = substr(line, 5, nchar(line)-4)

          if (nchar(unparsed_block_name)==8)
          {
            stop(paste("Invalid file: line ",line_counter,", new section of file needs a name: use /***name***/.",sep=""))
          }
          split_block_name = split_string_at_comma_ignoring_parentheses(unparsed_block_name)
          new_block_info = determine_block_type(split_block_name,blocks,line_counter,block_type,block_name,factor_number,transition_model_number,potential_function_number,importance_proposal_number,mh_proposal_number,independent_mh_proposal_number,m_proposal_number,enk_transform_number,transition_proposal_number,data_number,method_number)

          # expect input for each block in one of the following forms:
          # (a) /***evaluate_log_prior***/, followed by a C++ function
          # (b) /***evaluate_log_prior***/, in the case of a model with multiple factors, where this is the nth factor
          # (c) /***prior,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1
          # (d) /***prior,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1, in the case of a model with multiple factors, where this is the nth factor
          # (e) /***evaluate_log_prior,dlnorm(tau,1,1,TRUE)***/, so that we extract the named argument from the parameters read into the function, and input the other arguments, in this order, in to dlnorm (dlnorm in this case is a base R function, but it could be a user defined one in an R, python or julia file)
          # (f) /***evaluate_log_prior,dlnorm(tau,1,1,TRUE)***/, is the same case as (e), except that we have a model with multiple factors and this is the nth factor
          # (g) /***evaluate_log_likelihood,n,dnorm(data.y,0,recip(sqrt(tau)),TRUE)***/
          # (h) /***simulate_prior,tau=rlnorm(1,1,1)***/, when we simulate a variable
          # (i) /***simulate_prior,theta=rlnorm(1,1,1)***/, for simulations that involve multiple steps
          # (j) /***simulate_prior,tau=rlnorm(1,theta,1)***/, the counterpart to (i)

          block_type = new_block_info[[1]]
          number_to_pass_to_extract_block = new_block_info[[2]]
          block_name = new_block_info[[3]]
          is_custom = new_block_info[[4]]
          block_function = new_block_info[[5]]
          factor_number = new_block_info[[6]]
          transition_model_number = new_block_info[[7]]
          potential_function_number = new_block_info[[8]]
          importance_proposal_number = new_block_info[[9]]
          mh_proposal_number = new_block_info[[10]]
          independent_mh_proposal_number = new_block_info[[11]]
          m_proposal_number = new_block_info[[12]]
          enk_transform_number = new_block_info[[13]]
          transition_proposal_number = new_block_info[[14]]
          data_number = new_block_info[[15]]
          #reinforce_gradient_number = new_block_info[[17]]
          method_number = new_block_info[[16]]

        }
      }
    }

    if ( (in_block==TRUE) && (starting_block_flag==FALSE) && (is_custom==TRUE) )
    {
      block_code = paste(block_code,line,sep="\n")
    }

    # if (in_block==FALSE)
    # {
    #   header = paste(header,line,sep="\n")
    # }

  }

  if (in_block==TRUE)
  {
    # ignore block if block number is not positive
    if (number_to_pass_to_extract_block>0)
    {
      blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,parameter_list,external_packages,julia_bin_dir,julia_required_libraries)
      if ( (factor_number!=0) && (factor_number==length(blocks[["factor"]])) )
      {
        print_factor_info(length(blocks[["factor"]]),blocks,line_counter)
      }
      if ( (transition_model_number!=0) && (transition_model_number==length(blocks[["transition_model"]])) )
      {
        print_transition_model_info(length(blocks[["transition_model"]]),blocks,line_counter)
      }
      if ( (potential_function_number!=0) && (potential_function_number==length(blocks[["potential_function"]])) )
      {
        print_potential_function_info(length(blocks[["potential_function"]]),blocks,line_counter)
      }
      if ( (importance_proposal_number!=0) && importance_proposal_number==length(blocks[["importance_proposal"]]))
      {
        print_importance_proposal_info(length(blocks[["importance_proposal"]]),blocks,line_counter)
      }
      if ( (mh_proposal_number!=0) && mh_proposal_number==length(blocks[["mh_proposal"]]))
      {
        print_mh_proposal_info(length(blocks[["mh_proposal"]]),blocks,line_counter)
      }
      if ( (independent_mh_proposal_number!=0) && independent_mh_proposal_number==length(blocks[["independent_mh_proposal"]]))
      {
        print_independent_mh_proposal_info(length(blocks[["independent_mh_proposal"]]),blocks,line_counter)
      }
      if ( (m_proposal_number!=0) && m_proposal_number==length(blocks[["m_proposal"]]))
      {
        print_m_proposal_info(length(blocks[["m_proposal"]]),blocks,line_counter)
      }
      if ( (enk_transform_number!=0) && (enk_transform_number==length(blocks[["enk_transform"]])) )
      {
        print_enk_transform_info(length(blocks[["enk_transform"]]),blocks,line_counter)
      }
      if ( (transition_proposal_number!=0) && (transition_proposal_number==length(blocks[["transition_proposal"]])) )
      {
        print_transition_proposal_info(length(blocks[["transition_proposal"]]),blocks,line_counter)
      }
      if ( (data_number!=0) && (data_number==length(blocks[["data"]])) )
      {
        print_data_info(length(blocks[["data"]]),blocks,line_counter)
      }
      #if ( (reinforce_gradient_number!=0) && (reinforce_gradient_number==length(blocks[["reinforce_gradient"]])) )
      #{
      #  print_reinforce_gradient_info(length(blocks[["reinforce_gradient"]]),blocks,line_counter)
      #}
      if ( (method_number!=0) && (method_number==length(blocks[["method"]])) )
      {
        print_method_info(length(blocks[["method"]]),blocks,line_counter)
      }

    }
  }

  close(the_file)

  #browser()

  blocks = check_types(blocks)

  return(blocks)
}

# expect input for each block in one of the following forms:
# (a) /***evaluate_log_prior***/, followed by a C++ function
# (b) /***evaluate_log_prior,n***/, in the case of a model with multiple factors, where this is the nth factor
# (c) /***prior,norm,tau,0,1***/, to use a normal distribution with parameters 0 and 1
# (d) /***prior,n,norm,tau,0,1***/, to use a normal distribution with parameters 0 and 1, in the case of a model with multiple factors, where this is the nth factor

# expect input for each block in one of the following forms:
# (a) /***evaluate_log_prior***/, followed by a C++ function
# (b) /***evaluate_log_prior,n***/, in the case of a model with multiple factors, where this is the nth factor
# (c) /***prior,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1
# (d) /***prior,n,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1, in the case of a model with multiple factors, where this is the nth factor
# (e) /***evaluate_log_prior,dlnorm(tau,1,1,TRUE)***/, so that we extract the named argument from the parameters read into the function, and input the other arguments, in this order, in to dlnorm (dlnorm in this case is a base R function, but it could be a user defined one in an R, python or julia file)
# (f) /***evaluate_log_prior,n,dlnorm(tau,1,1,TRUE)***/, is the same case as (e), except that we have a model with multiple factors and this is the nth factor
# (g) /***evaluate_log_likelihood,n,dnorm(data$y,0,recip(sqrt(tau)),TRUE)***/
# (h) /***simulate_prior,tau=rlnorm(1,1,1)***/, when we simulate a variable
# (i) /***simulate_prior,1,theta=rlnorm(1,theta,1)***/, for simulations that involve multiple steps
# (j) /***simulate_prior,2,tau=rlnorm(1,theta,1)***/, the counterpart to (i)
