get_text_between_first_pair_of_matching_parentheses <- function(my_string)
{
  # Find the positions of open and closed brackets
  bracket_positions = c(0,0)
  bracket_positions[1] = regexpr("\\(", my_string)[1]
  bracket_positions[2] = regexpr("\\)", my_string)[1]

  # Check if there is at least one pair of brackets
  if (any(bracket_positions > 0))
  {
    text_between_parentheses = substr(my_string, bracket_positions[1], bracket_positions[2])
    return(text_between_parentheses)#return(gsub("\\(|\\)", "", text_between_parentheses))
  }
  else
  {
    return("")
  }
}

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
                        model_parameter_list = list())
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
        if (parameter_number<=length(model_parameter_list))
        {
          parameter_arguments = c(parameter_arguments,h)
          do.call("<-",list(h, model_parameter_list[[parameter_number]]))
        }
        else
        {
          stop(paste("In call ",input,", parameter ",as.numeric(substr(h,2,nchar(h)))," not found in model_parameter_list (model_parameter_list has length ",length(model_parameter_list),").",sep=""))
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
      else if ( ("data_variable" %in% names(current_factor_info)) && (block_name=="data_variable") )
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
      else if ( ("linear_gaussian_data_variable" %in% names(current_factor_info)) && (block_name=="linear_gaussian_data_variable") )
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("linear_gaussian_data_state_variable" %in% names(current_factor_info)) && (block_name=="linear_gaussian_data_state_variable") )
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
      else if ( ("nonlinear_gaussian_data_variable" %in% names(current_factor_info)) && (block_name=="nonlinear_gaussian_data_variable") )
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

  all_names = c(ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_types)

  # Get the current transition_model info.
  if ("transition_model" %in% names(blocks))
  {
    current_transition_model_info = blocks[["transition_model"]][[transition_model_number]]

    current_transition_model_names = names(current_transition_model_info)

    if (block_name %in% ilike_transition_model_types)
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
    else if (block_name %in% custom_transition_model_types)
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
      else if ( ("linear_gaussian_transition_variable" %in% names(current_transition_model_info)) && (block_name=="linear_gaussian_transition_variable") )
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
      else if ( ("nonlinear_gaussian_transition_variable" %in% names(current_transition_model_info)) && (block_name=="nonlinear_gaussian_transition_variable") )
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
        stop(paste("Invalid file: line ",line_counter,', mh_proposal specified, but previous mh proposal block is incomplete (did you specify both evaluate_log_mh_proposal and simulate_mh_proposal, and if proposing on a transformed space, the transform, inverse transform and jacobian_matrix?)',sep=""))
      }

    }
  }
  else
  {
    mh_proposal_number = mh_proposal_number + 1
  }

  return(mh_proposal_number)
}

unadjusted_proposal_processing = function(unadjusted_proposal_number,blocks,block_name,line_counter)
{
  # Is this a continuation of the current unadjusted_proposal, or a new one?

  # Get the current unadjusted_proposal info.
  if ("unadjusted_proposal" %in% names(blocks))
  {
    current_unadjusted_proposal_info = blocks[["unadjusted_proposal"]][[unadjusted_proposal_number]]

    if ("unadjusted_proposal" %in% names(current_unadjusted_proposal_info))
    {
      # unadjusted_proposal is complete
      print_unadjusted_proposal_info(unadjusted_proposal_number,blocks,line_counter-1)
      unadjusted_proposal_number = unadjusted_proposal_number + 1
    }
    else if ( ("simulate_unadjusted_proposal" %in% names(current_unadjusted_proposal_info)) || ("unadjusted_transform" %in% names(current_unadjusted_proposal_info)) || ("unadjusted_inverse_transform" %in% names(current_unadjusted_proposal_info)) || ("unadjusted_transform_jacobian_matrix" %in% names(current_unadjusted_proposal_info)) )
    {
      if ( ("unadjusted_transform" %in% names(current_unadjusted_proposal_info)) || ("unadjusted_inverse_transform" %in% names(current_unadjusted_proposal_info)) || ("unadjusted_transform_jacobian_matrix" %in% names(current_unadjusted_proposal_info)) )
      {
        # unadjusted_proposal is complete
        if ( ("simulate_unadjusted_proposal" %in% names(current_unadjusted_proposal_info)) && ("unadjusted_transform" %in% names(current_unadjusted_proposal_info)) && ("unadjusted_inverse_transform" %in% names(current_unadjusted_proposal_info)) && ("unadjusted_transform_jacobian_matrix" %in% names(current_unadjusted_proposal_info)) )
        {
          print_unadjusted_proposal_info(unadjusted_proposal_number,blocks,line_counter-1)
          unadjusted_proposal_number = unadjusted_proposal_number + 1
        }
      }
      else
      {
        # unadjusted_proposal is complete
        if ( ("simulate_unadjusted_proposal" %in% names(current_unadjusted_proposal_info)) )
        {
          print_unadjusted_proposal_info(unadjusted_proposal_number,blocks,line_counter-1)
          unadjusted_proposal_number = unadjusted_proposal_number + 1
        }
      }

      if (block_name=="unadjusted_proposal")
      {
        stop(paste("Invalid file: line ",line_counter,', unadjusted_proposal specified, but previous proposal block is incomplete (if proposing on a transformed space did you specify the transform, inverse transform and jacobian_matrix?)',sep=""))
      }

    }
  }
  else
  {
    unadjusted_proposal_number = unadjusted_proposal_number + 1
  }

  return(unadjusted_proposal_number)
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
        stop(paste("Invalid file: line ",line_counter,', transition_proposal specified, but previous transition proposal block is incomplete (did you specify both evaluate_log_transition_proposal and simulate_transition_proposal?)',sep=""))
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
    if (!grepl("XPtr",last_factor_names[j]))
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
    if (!grepl("XPtr",last_transition_model_names[j]))
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
    if (!grepl("XPtr",last_potential_function_names[j]))
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
    if (!grepl("XPtr",last_importance_proposal_names[j]))
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
    if (!grepl("XPtr",last_mh_proposal_names[j]))
      mh_proposal_info_string = paste(mh_proposal_info_string,last_mh_proposal_names[j],sep=", ")
  }
  if (nchar(mh_proposal_info_string)>2)
  {
    mh_proposal_info_string = substr(mh_proposal_info_string,3,nchar(mh_proposal_info_string))
  }
  print(paste('mh_proposal ends on line ',line_counter,'. Contains ',mh_proposal_info_string,'.',sep = ""))
}

print_unadjusted_proposal_info = function(unadjusted_proposal_index,blocks,line_counter)
{
  unadjusted_proposal_info_string = ""
  last_unadjusted_proposal_names = names(blocks[["unadjusted_proposal"]][[unadjusted_proposal_index]])
  for (j in 1:length(last_unadjusted_proposal_names))
  {
    if (!grepl("XPtr",last_unadjusted_proposal_names[j]))
      unadjusted_proposal_info_string = paste(unadjusted_proposal_info_string,last_unadjusted_proposal_names[j],sep=", ")
  }
  if (nchar(unadjusted_proposal_info_string)>2)
  {
    unadjusted_proposal_info_string = substr(unadjusted_proposal_info_string,3,nchar(unadjusted_proposal_info_string))
  }
  print(paste('unadjusted_proposal ends on line ',line_counter,'. Contains ',unadjusted_proposal_info_string,'.',sep = ""))
}

print_independent_mh_proposal_info = function(independent_mh_proposal_index,blocks,line_counter)
{
  independent_mh_proposal_info_string = ""
  last_independent_mh_proposal_names = names(blocks[["independent_mh_proposal"]][[independent_mh_proposal_index]])
  for (j in 1:length(last_independent_mh_proposal_names))
  {
    if (!grepl("XPtr",last_independent_mh_proposal_names[j]))
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
    if (!grepl("XPtr",last_m_proposal_names[j]))
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
    if (!grepl("XPtr",last_transition_proposal_names[j]))
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
    if (!grepl("XPtr",last_enk_transform_names[j]))
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
    if (!grepl("XPtr",last_reinforce_gradient_names[j]))
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
    if (!grepl("XPtr",last_method_names[j]))
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
    if (!grepl("XPtr",last_data_names[j]))
      data_info_string = paste(data_info_string,last_data_names[j],sep=", ")
  }
  if (nchar(data_info_string)>2)
  {
    data_info_string = substr(data_info_string,3,nchar(data_info_string))
  }
  print(paste('data ends on line ',line_counter,'. Contains ',data_info_string,'.',sep = ""))
}

determine_block_type = function(split_block_name,blocks,line_counter,block_type,block_name,factor_number,transition_model_number,potential_function_number,importance_proposal_number,mh_proposal_number,unadjusted_proposal_number,independent_mh_proposal_number,m_proposal_number,enk_transform_number,transition_proposal_number,data_number,method_number)
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
  sbi_likelihood_function_types = c("simulate_data_model","sbi_likelihood","summary_statistics","data_variable")
  linear_gaussian_data_model_types = c("linear_gaussian_data_model","linear_gaussian_data_matrix","linear_gaussian_data_covariance","linear_gaussian_data_variable","linear_gaussian_data_state_variable")
  nonlinear_gaussian_data_model_types = c("nonlinear_gaussian_data_model","nonlinear_gaussian_data_function","nonlinear_gaussian_data_covariance","nonlinear_gaussian_data_variable")
  other_likelihood_function_types = c("algorithmic_likelihood")
  factor_function_types = c(prior_function_types,custom_likelihood_function_types,sbi_likelihood_function_types,linear_gaussian_data_model_types,nonlinear_gaussian_data_model_types,other_likelihood_function_types)
  ilike_transition_model_types = c("transition_model")
  linear_gaussian_transition_model_types = c("linear_gaussian_transition_model","linear_gaussian_transition_matrix","linear_gaussian_transition_covariance","linear_gaussian_transition_variable")
  nonlinear_gaussian_transition_model_types = c("nonlinear_gaussian_transition_model","nonlinear_gaussian_transition_function","nonlinear_gaussian_transition_covariance","nonlinear_gaussian_transition_variable")
  custom_transition_model_types = c("simulate_transition_model","evaluate_log_transition_model")
  transition_model_types = c(ilike_transition_model_types,linear_gaussian_transition_model_types,nonlinear_gaussian_transition_model_types,custom_transition_model_types)
  custom_potential_function_types = c("evaluate_log_potential_function")
  ilike_potential_function_types = c("potential_function")
  potential_function_types = c(custom_potential_function_types,ilike_potential_function_types)
  data_function_types = c("data")
  importance_proposal_types = c("simulate_importance_proposal","evaluate_log_importance_proposal","importance_proposal")
  mh_proposal_types = c("simulate_mh_proposal","evaluate_log_mh_proposal","mh_proposal","mh_transform","mh_inverse_transform","mh_transform_jacobian_matrix","mh_factor_index")
  unadjusted_proposal_types = c("simulate_unadjusted_proposal","unadjusted_proposal","unadjusted_transform","unadjusted_inverse_transform","unadjusted_transform_jacobian_matrix","unadjusted_factor_index")
  independent_mh_proposal_types = c("simulate_independent_mh_proposal","evaluate_log_independent_mh_proposal","independent_mh_proposal","independent_mh_transform","independent_mh_inverse_transform","independent_mh_transform_jacobian_matrix","independent_mh_factor_index")
  m_proposal_types = c("simulate_m_proposal","m_proposal","m_transform","m_inverse_transform","m_transform_jacobian_matrix","m_factor_index")
  transition_proposal_types = c("simulate_transition_proposal","evaluate_log_transition_proposal")
  enk_transform_types = c("enk_transform","enk_inverse_transform")
  #reinforce_gradient = c("reinforce_gradient")
  method_function_types = c("mcmc_weights","mcmc_termination","adaptive_resampling","adaptive_target","smc_termination","smc_sequence","reinforce_gradient","enk_likelihood_index","enk_shifter","filter")

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
  else if (block_name %in% unadjusted_proposal_types)
  {
    block_type = "unadjusted_proposal"
    unadjusted_proposal_number = unadjusted_proposal_processing(unadjusted_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = unadjusted_proposal_number
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
              unadjusted_proposal_number,
              independent_mh_proposal_number,
              m_proposal_number,
              enk_transform_number,
              transition_proposal_number,
              data_number,
              method_number))
}

get_parameters_output_function_body <- function(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
{
  if (output_variable=="")
  {
    #stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))

    if (R_functions==TRUE)
    {
      function_body = paste('function(',R_args,') { return(',R_function_name,'(',R_function_arguments_string,'))}',sep="")
    }
    else
    {
      function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); return list_to_parameters(f(',cpp_function_arguments_string,'));',sep="")
    }
  }
  else
  {
    if (R_functions==TRUE)
    {
      function_body = paste('function(',R_args,') { output=list(); output[["',output_variable,'"]] = ',R_function_name,'(',R_function_arguments_string,'); return(output)}',sep="")
    }
    else
    {
      function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = Rcpp::as<arma::mat>(NumericMatrix(f(',cpp_function_arguments_string,'))); return output;',sep="")
    }
  }

  return(function_body)
}

get_double_output_function_body <- function(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
{
  if (R_functions==TRUE)
  {
    function_body = paste('function(',R_args,') { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
  }
  else
  {
    function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
  }

  return(function_body)
}

get_matrix_output_function_body <- function(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
{
  if (R_functions==TRUE)
  {
    function_body = paste('function(',R_args,') { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
  }
  else
  {
    function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','NumericMatrix num_mat = NumericMatrix(f(',cpp_function_arguments_string,')); return arma::mat(num_mat.begin(), num_mat.nrow(), num_mat.ncol(), false);',sep="")
  }

  return(function_body)
}

get_list_output_function_body <- function(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
{
  if (R_functions==TRUE)
  {
    function_body = paste('function(',R_args,') { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
  }
  else
  {
    function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','List list = List(f(',cpp_function_arguments_string,')); return list;',sep="")
  }

  return(function_body)
}

get_string_output_function_body <- function(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
{
  if (R_functions==TRUE)
  {
    function_body = paste('function(',R_args,') { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
  }
  else
  {
    function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return Rcpp::as<std::string>(f(',cpp_function_arguments_string,'));',sep="")
  }

  return(function_body)
}

parameter_types_for_cpp <- function(block_code)
{
  #arguments = regmatches(block_code, gregexpr("\\([^)]+\\)", block_code))[[1]][1]
  #arguments = gsub("\\(|\\)", "", arguments)
  arguments = get_text_between_first_pair_of_matching_parentheses(block_code)
  argument_vector = unlist(strsplit(arguments, ","))

  if (!is.null(argument_vector))
  {
    which_proposed_parameters = matrix(0,length(argument_vector))
    which_parameters = matrix(0,length(argument_vector))
    which_proposal_parameters = matrix(0,length(argument_vector))
    which_data = matrix(0,length(argument_vector))

    for (i in 1:length(argument_vector))
    {

      if (grepl("&", argument_vector[i]))
      {
        after_ampersand = unlist(strsplit(argument_vector[i], "&"))
        parameter_name = trimws(after_ampersand[length(after_ampersand)])
        parameter_name = gsub("\\(|\\)", "", parameter_name)

        if (parameter_name=="proposed_parameters")
        {
          which_proposed_parameters[i] = 1
        }

        if (parameter_name=="parameters")
        {
          which_parameters[i] = 1
        }

        if (parameter_name=="proposal_parameters")
        {
          which_proposal_parameters[i] = 1
        }

        if (parameter_name=="data")
        {
          which_data[i] = 1
        }
      }
    }
  }
  else
  {
    which_proposed_parameters = NULL
    which_parameters = NULL
    which_proposal_parameters = NULL
    which_data = NULL
  }

  return(list(which_proposed_parameters,which_parameters,which_proposal_parameters,which_data))
}

return_types_for_cpp <- function(block_name,line_counter,which_proposed_parameters,which_parameters,which_proposal_parameters,which_data)
{
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
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a data function.",sep=""))
    }

    return_type = "Data"
    arguments = c()
    #args_for_typedef = ""
    R_args = ""
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
    R_args = "parameters"
  }
  else if (block_name=="evaluate_gradient_log_prior")
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
    arguments[1] = "const std::string &variable"
    arguments[2] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "variable,parameters"
  }
  else if (block_name=="evaluate_second_gradient_log_prior")
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
    arguments[1] = "const std::string &variable1"
    arguments[2] = "const std::string &variable2"
    arguments[3] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"

    R_args = "variable1,variable2,parameters"
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

    R_args = "parameters,data"
  }
  else if (block_name=="evaluate_gradient_log_likelihood")
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
    arguments[1] = "const std::string &variable"
    arguments[2] = "const Parameters &parameters"
    arguments[3] = "const Data &data"
    #args_for_typedef = "const Parameters&"

    R_args = "variable,parameters,data"
  }
  else if (block_name=="evaluate_second_gradient_log_likelihood")
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
    arguments[1] = "const std::string &variable1"
    arguments[2] = "const std::string &variable2"
    arguments[3] = "const Parameters &parameters"
    arguments[4] = "const Data &data"
    #args_for_typedef = "const Parameters&"

    R_args = "variable1,variable2,parameters,data"
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "RandomNumberGenerator &rng"
    #args_for_typedef = "RandomNumberGenerator&"

    R_args = ""

  }
  else if (block_name=="evaluate_log_importance_proposal")
  {

    if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0)  )
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Parameters &proposal_parameters"
      arguments[3] = "const Data &data"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,proposal_parameters,data"
      proposal_type = 4
    }
    else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Data &data"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,data"
      proposal_type = 3
    }
    else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Parameters &proposal_parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,proposal_parameters"
      proposal_type = 2
    }
    else
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters"
      proposal_type = 1
    }
  }
  else if (block_name=="simulate_importance_proposal")
  {

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
      R_args = "proposal_parameters,data"

      proposal_type = 4
    }
    else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
    {
      arguments = c()
      arguments[1] = "RandomNumberGenerator &rng"
      arguments[2] = "const Data &data"
      #args_for_typedef = "const Parameters&"
      R_args = "data"

      proposal_type = 3
    }
    else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
    {
      arguments = c()
      arguments[1] = "RandomNumberGenerator &rng"
      arguments[2] = "const Parameters &proposal_parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "proposal_parameters,data"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      arguments[1] = "RandomNumberGenerator &rng"
      #args_for_typedef = "const Parameters&"
      R_args = ""
      proposal_type = 1
    }
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
      R_args = "proposed_parameters,parameters,proposal_parameters,data"

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
      R_args = "proposed_parameters,parameters,data"
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
      R_args = "proposed_parameters,parameters,proposal_parameters"
      proposal_type = 2
    }
    else
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,parameters"
      proposal_type = 1
    }
  }
  else if (block_name=="simulate_mh_proposal")
  {
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
        R_args = "parameters,proposal_parameters,data"
        proposal_type = 4
      }
      else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Data &data"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,data"
        proposal_type = 3
      }
      else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Parameters &proposal_parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,proposal_parameters"
        proposal_type = 2
      }
      else
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters"
        proposal_type = 1
      }
    }
    else
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, should this have been specified as simulate_independent_mh_proposal?",sep=""))
    }
  }
  else if (block_name=="simulate_unadjusted_proposal")
  {
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
        R_args = "parameters,proposal_parameters,data"
        proposal_type = 4
      }
      else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Data &data"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,data"
        proposal_type = 3
      }
      else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Parameters &proposal_parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,proposal_parameters"
        proposal_type = 2
      }
      else
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters"
        proposal_type = 1
      }
    }
    else
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, should this have been specified as simulate_independent_mh_proposal?",sep=""))
    }
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
      R_args = "proposed_parameters,proposal_parameters,data"
      proposal_type = 4
    }
    else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Data &data"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,data"
      proposal_type = 3
    }
    else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      arguments[2] = "const Parameters &proposal_parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "proposed_parameters,proposal_parameters"
      proposal_type = 2
    }
    else
    {
      return_type = "double"
      arguments = c()
      arguments[1] = "const Parameters &proposed_parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 1
    }
  }
  else if (block_name=="simulate_independent_mh_proposal")
  {
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
        R_args = "proposal_parameters,data"
        proposal_type = 4
      }
      else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Data &data"
        #args_for_typedef = "const Parameters&"
        R_args = "data"
        proposal_type = 3
      }
      else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &proposal_parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "proposal_parameters"
        proposal_type = 2
      }
      else
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        #args_for_typedef = "const Parameters&"
        R_args = ""
        proposal_type = 1
      }
    }
    else
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": parameters used, should this have been specified as simulate_mh_proposal?",sep=""))
    }
  }
  else if (block_name=="simulate_m_proposal")
  {
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
        R_args = "parameters,proposal_parameters,data"
        proposal_type = 4
      }
      else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Data &data"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,data"
        proposal_type = 3
      }
      else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        arguments[3] = "const Parameters &proposal_parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters,proposal_parameters"
        proposal_type = 2
      }
      else
      {
        arguments = c()
        arguments[1] = "RandomNumberGenerator &rng"
        arguments[2] = "const Parameters &parameters"
        #args_for_typedef = "const Parameters&"
        R_args = "parameters"
        proposal_type = 1
      }
    }
    else
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, invalid proposal type?",sep=""))
    }

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
  }
  else if (block_name=="enk_likelihood_index")
  {
    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a enk_likelihood_index function.",sep=""))
    }

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a enk_likelihood_index function.",sep=""))
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
  }
  else if ( (block_name=="m_factor_index") || (block_name=="mh_factor_index") || (block_name=="independent_mh_factor_index") || (block_name=="unadjusted_factor_index") )
  {
    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a factor_index function.",sep=""))
    }

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a factor_index function.",sep=""))
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
  }
  else if (block_name=="simulate_data_model")
  {
    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a data_model.",sep=""))
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
    R_args = "parameters"
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

    return_type = "Data"
    arguments = c()
    arguments[1] = "const Data &data"
    #args_for_typedef = "const Parameters&"
    R_args = "data"

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

    return_type = "Data"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "RandomNumberGenerator &rng"
    arguments[2] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "RandomNumberGenerator &rng"
    arguments[2] = "const Parameters &parameters"
    arguments[3] = "const Data &data"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters,data"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
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

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"

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

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"
  }
  else if ( (block_name=="linear_gaussian_data_variable") || (block_name=="linear_gaussian_data_state_variable") )
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

    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
    }

    arguments = c()
    R_args = ""

    return_type = "std::string"

  }
  else if (block_name=="data_variable")
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

    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
    }

    arguments = c()
    R_args = ""

    return_type = "std::string"

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

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"
  }
  else if (block_name=="nonlinear_gaussian_data_variable")
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

    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
    }

    arguments = c()
    R_args = ""

    return_type = "std::string"

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

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"
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

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"

  }
  else if (block_name=="linear_gaussian_transition_variable")
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

    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
    }

    arguments = c()
    R_args = ""

    return_type = "std::string"

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

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    if (sum(which_parameters)>0)
    {
      arguments = c()
      arguments[1] = "const Parameters &parameters"
      #args_for_typedef = "const Parameters&"
      R_args = "parameters"
      proposal_type = 2
    }
    else
    {
      arguments = c()
      R_args = ""
      proposal_type = 1
    }

    return_type = "arma::mat"
  }
  else if (block_name=="nonlinear_gaussian_transition_variable")
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

    if (sum(which_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
    }

    arguments = c()
    R_args = ""

    return_type = "std::string"

  }
  else if (block_name=="evaluate_log_transition_model")
  {

    if (sum(which_proposal_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
    }

    if (sum(which_data)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
    }

    return_type = "double"
    arguments = c()
    arguments[1] = "const Parameters &proposed_parameters"
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "proposed_parameters,parameters"
  }
  else if (block_name=="evaluate_log_potential_function")
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
    R_args = "parameters,data"
  }
  else if (block_name=="evaluate_log_transition_proposal")
  {

    if (sum(which_proposal_parameters)>0)
    {
      stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
    }

    return_type = "double"
    arguments = c()
    arguments[1] = "const Parameters &proposed_parameters"
    arguments[2] = "const Parameters &parameters"
    arguments[3] = "const Data &data"
    #args_for_typedef = "const Parameters&"
    R_args = "proposed_parameters,parameters,data"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
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
    R_args = "parameters"
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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
    R_args = "parameters"

  }
  else if (block_name=="unadjusted_transform")
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
  }
  else if (block_name=="unadjusted_inverse_transform")
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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"
  }
  else if (block_name=="unadjusted_transform_jacobian_matrix")
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
    R_args = "parameters"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

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

    return_type = "Parameters"
    arguments = c()
    arguments[1] = "const Parameters &parameters"
    #args_for_typedef = "const Parameters&"
    R_args = "parameters"

  }
  else if (block_name=="independent_mh_transform_jacobian_matrix")
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
    R_args = "parameters"

  }
  else
  {
    stop(paste('Block ',block_name,', line number ',line_counter,': block is of unknown type. Did you specify an ilike function? If so, include "ilike::" before the function name.',sep=""))
  }

  args_for_typedef = paste(arguments,collapse=",")

  return(list(return_type,args_for_typedef,R_args,proposal_type))
}

function_name_for_cpp <- function(block_code)
{
  extracted_text = sub("\\(.*", "", block_code)
  split_at_space = strsplit(extracted_text," ")[[1]]
  return(split_at_space[length(split_at_space)])
}

extract_block <- function(blocks,block_type,block_name,factor_number,line_counter,block_code,block_function,is_custom,model_parameter_list,R_functions,external_packages,julia_bin_dir,julia_required_libraries,fileConn,verify_cpp_function_types,nesting_level,keep_temporary_model_code)
{

  # Get information about the order in which MCMC moves are included.
  if ( (block_type=="mh_proposal") || (block_type=="unadjusted_proposal") || (block_type=="independent_mh_proposal") || (block_type=="m_proposal") )
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

      if (block_type=="unadjusted_proposal")
      {
        blocks[["order_of_mcmc"]] = c(blocks[["order_of_mcmc"]],4)
      }
    }
  }

  cpp_function_name = paste(block_name,factor_number,sep="")
  R_function_name = paste(cpp_function_name,'_R',sep="")

  #temp_filename = paste(cpp_function_name,".cpp",sep="")
  xptr_name = paste(cpp_function_name,"getXPtr",sep="_")

  # ignore block if block number is not positive
  if (factor_number>0)
  {
    if (is_custom==TRUE)
    {
      pt = parameter_types_for_cpp(block_code)
      which_proposed_parameters = pt[[1]]
      which_parameters = pt[[2]]
      which_proposal_parameters = pt[[3]]
      which_data = pt[[4]]

      actual_cpp_function_name = function_name_for_cpp(block_code)

      rt = return_types_for_cpp(block_name,line_counter,which_proposed_parameters,which_parameters,which_proposal_parameters,which_data)
      return_type = rt[[1]]
      args_for_typedef = rt[[2]]
      R_args = rt[[3]]
      proposal_type = rt[[4]]

      writeLines(c(paste("SEXP ",xptr_name,"();",sep=""),
                   '\n',
                   block_code,
                   '\n',
                   "// [[Rcpp::export]]",
                   paste("SEXP ",xptr_name,"() {",sep=""),
                   paste("  typedef", return_type, "(*funcPtr)(", args_for_typedef, ");"),
                   paste("  return XPtr<funcPtr>(new funcPtr(&", actual_cpp_function_name, "));"),
                   "}",
                   '\n'),
                 fileConn)

      if (verify_cpp_function_types)
      {
        my_list = list(xptr_name,0,RcppXPtrUtils::cppXPtr(block_code,plugins=c("cpp11"),depends = paste(c("ilike","RcppArmadillo","BH","dqrng","sitmo"),external_packages) ))
        names(my_list) = c(xptr_name,block_name,paste("verify_",xptr_name,sep=""))
      }
      else
      {
        my_list = list(xptr_name,0)
        names(my_list) = c(xptr_name,block_name)
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
        blocks[[block_type]][[factor_number]] = append(blocks[[block_type]][[factor_number]],my_list)
      }

      if (proposal_type!=0)
      {
        blocks[[block_type]][[factor_number]][["type"]] = proposal_type
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

          if (grepl("=", ilike_split_name[1]))
          {
            stop("When specifying an ilike function, you should not have an output (e.g. don't use x=ilike::mvnorm(...)).")
          }
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
        if ( (block_name=="prior") || (block_name=="importance_proposal") || (block_name=="independent_mh_proposal") || (block_name=="mh_proposal") || (block_name=="unadjusted_proposal") || (block_name=="m_proposal") || (block_name=="smc_sequence") || (block_name=="linear_gaussian_transition_model") || (block_name=="nonlinear_gaussian_transition_model") || (block_name=="linear_gaussian_data_model") || (block_name=="nonlinear_gaussian_data_model") )
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
        else if ( (block_name=="mcmc_termination") || (block_name=="adaptive_resampling") || (block_name=="adaptive_target") || (block_name=="smc_termination") || (block_name=="reinforce_gradient") || (block_name=="enk_likelihood_index") || (block_name=="enk_shifter") || (block_name=="filter") || (block_name=="sbi_likelihood") )
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
        else if (block_name=="algorithmic_likelihood")
        {
          split_arg_string = split_string_at_comma_ignoring_parentheses(arg_string)

          sub_filename = split_arg_string[1]

          parameters = list()
          if (length(split_arg_string)>1)
          {
            for (k in 2:(length(split_arg_string)))
            {
              parameters[[k-1]] = split_arg_string[k]
            }
          }

          sub_nesting_level = nesting_level + 1

          sub_filename = gsub('"', '', sub_filename)
          sub_filename = gsub("'", '', sub_filename)
          my_list = list(list(type=ilike_type,
                              model=ilike::compile(filename=sub_filename,model_parameter_list=model_parameter_list,external_packages=external_packages,julia_bin_dir=julia_bin_dir,julia_required_libraries=julia_required_libraries,verify_cpp_function_types=verify_cpp_function_types,keep_temporary_model_code=keep_temporary_model_code,nesting_level=sub_nesting_level),
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
                                    model_parameter_list)

        R_function_arguments = function_info[[2]]

        cpp_function_arguments_string = ""
        R_function_arguments_string = ""
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
            split_at_dollar = strsplit(R_function_arguments[i],"\\$")[[1]]

            any_lists = FALSE

            if (split_at_dollar[1]=="proposed_parameters")
            {
              which_proposed_parameters[i] = 1
              any_lists = TRUE
            }

            if (split_at_dollar[1]=="parameters")
            {
              which_parameters[i] = 1
              any_lists = TRUE
            }

            if (split_at_dollar[1]=="proposal_parameters")
            {
              which_proposal_parameters[i] = 1
              any_lists = TRUE
            }

            if (split_at_dollar[1]=="data")
            {
              which_data[i] = 1
              any_lists = TRUE
            }

            if (length(split_at_dollar)==0)
            {
              stop(paste("Block ",block_name,", line number ",line_counter,": argument is of size zero.",sep=""))
            }
            else if (length(split_at_dollar)==1)
            {
              if (any_lists)
              {
                cpp_function_arguments_string = paste(cpp_function_arguments_string,',parameters_to_list(',split_at_dollar[1],')',sep="")
              }
              else
              {
                cpp_function_arguments_string = paste(cpp_function_arguments_string,split_at_dollar[1],sep=",")
              }
              R_function_arguments_string = paste(R_function_arguments_string,split_at_dollar[1],sep=",")
            }
            else if (length(split_at_dollar)>=2)
            {
              cpp_function_arguments_string = paste(cpp_function_arguments_string,',',split_at_dollar[1],'["',paste(split_at_dollar[2:length(split_at_dollar)],collapse=""),'"]',sep="")
              R_function_arguments_string = paste(R_function_arguments_string,',',split_at_dollar[1],'$',paste(split_at_dollar[2:length(split_at_dollar)],collapse=""),sep="")
            }
          }

          cpp_function_arguments_string = substr(cpp_function_arguments_string,2,nchar(cpp_function_arguments_string))
          R_function_arguments_string = substr(R_function_arguments_string,2,nchar(R_function_arguments_string))
        }

        # give C++ wrapper for this
        # write to cpp file
        # source
        # delete file

        proposal_type = 0

        # if (R_functions)
        # {
        #   R_function_arguments_string = ""
        #   for (r in 1:length(R_function_arguments))
        #   {
        #     paste(R_function_arguments_string,sep="")
        #   }
        # }

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
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a data function.",sep=""))
          }

          return_type = "Data"
          arguments = c()
          #args_for_typedef = ""
          R_args = ""

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
          R_args = "parameters"

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="evaluate_gradient_log_prior")
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
          arguments[1] = "const std::string &variable"
          arguments[2] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "variable,parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="evaluate_second_gradient_log_prior")
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
          arguments[1] = "const std::string &variable1"
          arguments[2] = "const std::string &variable2"
          arguments[3] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"

          R_args = "variable1,variable2,parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          R_args = "parameters,data"

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="evaluate_gradient_log_likelihood")
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
          arguments[1] = "const std::string &variable"
          arguments[2] = "const Parameters &parameters"
          arguments[3] = "const Data &data"
          #args_for_typedef = "const Parameters&"

          R_args = "variable,parameters,data"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="evaluate_second_gradient_log_likelihood")
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
          arguments[1] = "const std::string &variable1"
          arguments[2] = "const std::string &variable2"
          arguments[3] = "const Parameters &parameters"
          arguments[4] = "const Data &data"
          #args_for_typedef = "const Parameters&"

          R_args = "variable1,variable2,parameters,data"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          #args_for_typedef = "RandomNumberGenerator&"

          R_args = ""

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
            R_args = "proposed_parameters,proposal_parameters,data"
            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters,data"
            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters,proposal_parameters"
            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters"
            proposal_type = 1
          }

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="simulate_importance_proposal")
        {
          # if (output_variable=="")
          # {
          #   stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
          # }

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
            R_args = "proposal_parameters,data"

            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"
            R_args = "data"

            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "proposal_parameters,data"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            #args_for_typedef = "const Parameters&"
            R_args = ""
            proposal_type = 1
          }

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
            R_args = "proposed_parameters,parameters,proposal_parameters,data"

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
            R_args = "proposed_parameters,parameters,data"
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
            R_args = "proposed_parameters,parameters,proposal_parameters"
            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters,parameters"
            proposal_type = 1
          }

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="simulate_mh_proposal")
        {

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
              R_args = "parameters,proposal_parameters,data"
              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,data"
              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,proposal_parameters"
              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters"
              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, should this have been specified as simulate_independent_mh_proposal?",sep=""))
          }

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="simulate_unadjusted_proposal")
        {

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
              R_args = "parameters,proposal_parameters,data"
              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,data"
              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,proposal_parameters"
              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters"
              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, should this have been specified as simulate_independent_mh_proposal?",sep=""))
          }

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
            R_args = "proposed_parameters,proposal_parameters,data"
            proposal_type = 4
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters,data"
            proposal_type = 3
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "proposed_parameters,proposal_parameters"
            proposal_type = 2
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &proposed_parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 1
          }

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="simulate_independent_mh_proposal")
        {
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
              R_args = "proposal_parameters,data"
              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Data &data"
              #args_for_typedef = "const Parameters&"
              R_args = "data"
              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "proposal_parameters"
              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              #args_for_typedef = "const Parameters&"
              R_args = ""
              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": parameters used, should this have been specified as simulate_mh_proposal?",sep=""))
          }

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="simulate_m_proposal")
        {

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
              R_args = "parameters,proposal_parameters,data"
              proposal_type = 4
            }
            else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Data &data"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,data"
              proposal_type = 3
            }
            else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              arguments[3] = "const Parameters &proposal_parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters,proposal_parameters"
              proposal_type = 2
            }
            else
            {
              arguments = c()
              arguments[1] = "RandomNumberGenerator &rng"
              arguments[2] = "const Parameters &parameters"
              #args_for_typedef = "const Parameters&"
              R_args = "parameters"
              proposal_type = 1
            }
          }
          else
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no parameters used, invalid proposal type?",sep=""))
          }

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          if (R_functions==TRUE)
          {
            function_body = paste('function() { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
          }
          else
          {
            function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); NumericVector output = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
          }

        }
        else if (block_name=="enk_likelihood_index")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a enk_likelihood_index function.",sep=""))
          }

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a enk_likelihood_index function.",sep=""))
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

          if (R_functions==TRUE)
          {
            function_body = paste('function() { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
          }
          else
          {
            function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); NumericVector output = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
          }

        }
        else if ( (block_name=="m_factor_index") || (block_name=="mh_factor_index") || (block_name=="independent_mh_factor_index") || (block_name=="unadjusted_factor_index") )
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a factor_index function.",sep=""))
          }

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a factor_index function.",sep=""))
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

          if (R_functions==TRUE)
          {
            function_body = paste('function() { return(',R_function_name,'(',R_function_arguments_string,')) }',sep="")
          }
          else
          {
            function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); NumericVector output = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
          }

        }
        else if (block_name=="simulate_data_model")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in a data_model.",sep=""))
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
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Data"
          arguments = c()
          arguments[1] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          R_args = "data"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Data"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          arguments[2] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          arguments[2] = "const Parameters &parameters"
          arguments[3] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters,data"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          if (sum(which_parameters)>0)
          {
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            R_args = ""
            proposal_type = 1
          }

          return_type = "arma::mat"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
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

          if (sum(which_parameters)>0)
          {
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            R_args = ""
            proposal_type = 1
          }

          return_type = "arma::mat"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if ( (block_name=="linear_gaussian_data_variable") || (block_name=="linear_gaussian_data_state_variable") )
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

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          arguments = c()
          R_args = ""

          return_type = "std::string"

          function_body = get_string_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
        }
        else if (block_name=="data_variable")
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

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          arguments = c()
          R_args = ""

          return_type = "std::string"

          function_body = get_string_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
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
          R_args = "parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="nonlinear_gaussian_data_variable")
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

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          arguments = c()
          R_args = ""

          return_type = "std::string"

          function_body = get_string_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
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

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_parameters)>0)
          {
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            R_args = ""
            proposal_type = 1
          }

          return_type = "arma::mat"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_parameters)>0)
          {
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            R_args = ""
            proposal_type = 1
          }

          return_type = "arma::mat"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="linear_gaussian_transition_variable")
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

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          arguments = c()
          R_args = ""

          return_type = "std::string"

          function_body = get_string_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
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

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          if (sum(which_parameters)>0)
          {
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
            R_args = "parameters"
            proposal_type = 2
          }
          else
          {
            arguments = c()
            R_args = ""
            proposal_type = 1
          }

          return_type = "arma::mat"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="nonlinear_gaussian_transition_variable")
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

          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in this function.",sep=""))
          }

          arguments = c()
          R_args = ""

          return_type = "std::string"

          function_body = get_string_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
        }
        else if (block_name=="evaluate_log_transition_model")
        {

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from data not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &proposed_parameters"
          arguments[2] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "proposed_parameters,parameters"

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
          R_args = "parameters,data"

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="evaluate_log_transition_proposal")
        {

          if (sum(which_proposal_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from proposal parameters not possible in this function.",sep=""))
          }

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &proposed_parameters"
          arguments[2] = "const Parameters &parameters"
          arguments[3] = "const Data &data"
          #args_for_typedef = "const Parameters&"
          R_args = "proposed_parameters,parameters,data"

          function_body = get_double_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)
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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
          R_args = "parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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
          R_args = "parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="unadjusted_transform")
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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="unadjusted_inverse_transform")
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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="unadjusted_transform_jacobian_matrix")
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
          R_args = "parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

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

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          R_args = "parameters"

          function_body = get_parameters_output_function_body(output_variable,R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else if (block_name=="independent_mh_transform_jacobian_matrix")
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
          R_args = "parameters"

          function_body = get_matrix_output_function_body(R_functions,R_args,R_function_name,R_function_arguments_string,cpp_function_arguments_string)

        }
        else
        {
          stop(paste('Block ',block_name,', line number ',line_counter,': block is of unknown type. Did you specify an ilike function? If so, include "ilike::" before the function name.',sep=""))
        }

        if (R_functions)
        {
          eval(parse(text=paste(R_function_name,function_info[[1]],sep="")))
          my_list = list(eval(parse(text=function_body)))
          names(my_list) = c(block_name)
        }
        else
        {
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

          writeLines(c('/*** R',
                       paste(R_function_name,function_info[[1]],sep=""),
            '*/',
            '\n',
            'using namespace Rcpp;',
            '#include <ilike.h>',
            '\n',
            paste("SEXP ",xptr_name,"();",sep=""),
            '\n',
            code,
            '\n',
            "// [[Rcpp::export]]",
            paste("SEXP ",xptr_name,"() {",sep=""),
            paste("  typedef", return_type, "(*funcPtr)(", args_for_typedef, ");"),
            paste("  return XPtr<funcPtr>(new funcPtr(&", cpp_function_name, "));"),
            "}",
            '\n'),
            fileConn)

          # writeLines(c(
          #   '/*** R',
          #   paste(R_function_name,function_info[[1]],sep=""),
          #   '*/'), fileConn)

          #my_list = list(get(xptr_name)())
          my_list = list(xptr_name,0)
          names(my_list) = c(xptr_name,block_name)
        }
        #my_list = list(RcppXPtrUtils::cppXPtr(code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo")))

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

        if (proposal_type!=0)
        {
          my_list[["type"]] = proposal_type

          if ("type" %in% names(blocks[[block_type]][[factor_number]]))
          {
            if (blocks[[block_type]][[factor_number]][["type"]]!=proposal_type)
              stop(paste("Block ",block_name,", line number ",line_counter,": different parts of the proposal/transition model/measurement model are of different types (maybe one relies on additional parameters and/or data, and the other does not?)",sep=""))
          }
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
      current_factor = blocks[["data"]][[i]]
      if ("data" %in% names(current_factor))
      {
        which_index = which(names(current_factor) == "data")
        for (j in 1:length(which_index))
        {
          if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
            RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Data")
        }
      }
    }
  }

  for (i in 1:length(blocks[["factor"]]))
  {
    current_factor = blocks[["factor"]][[i]]

    if ("evaluate_log_prior" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_prior")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "double", c("const Parameters&"))
      }
      }
    }

    if ("simulate_prior" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_prior")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&"))
      }
      }
    }

    if ("evaluate_log_likelihood" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_likelihood")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Data&"))
      }
      }
    }

    if ("simulate_data_model" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_data_model")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "Data", c("RandomNumberGenerator&","const Parameters&"))
      }
      }
    }

    if ("summary_statistics" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "summary_statistics")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "Data", c("const Data&"))
      }
      }
    }

    if ("nonlinear_gaussian_data_function" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "nonlinear_gaussian_data_function")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "Data", c("const Parameters&"))
      }
      }
    }

    if ("linear_gaussian_data_matrix" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_data_matrix")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c()) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid linear_gaussian_data_matrix specified.")
        }

        if ("type" %in% names(blocks[["factor"]][[i]]))
        {
          if (blocks[["factor"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("linear_gaussian_data_matrix and linear_gaussian_data_covariance functions are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both linear_gaussian_data_matrix and linear_gaussian_data_covariance.")
        }
        #blocks[["linear_gaussian_data_matrix"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("linear_gaussian_data_covariance" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_data_covariance")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c()) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid linear_gaussian_data_covariance specified.")
        }

        if ("type" %in% names(blocks[["factor"]][[i]]))
        {
          if (blocks[["factor"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("linear_gaussian_data_matrix and linear_gaussian_data_covariance functions are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both linear_gaussian_data_matrix and linear_gaussian_data_covariance.")
        }

        #blocks[["linear_gaussian_data_covariance"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("linear_gaussian_data_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_data_variable")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
      }
      }
    }

    if ("linear_gaussian_data_state_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_data_variable")
      for (j in 1:length(which_index))
      {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
        }
      }
    }

    if ("nonlinear_gaussian_data_covariance" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "nonlinear_gaussian_data_covariance")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c()) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid nonlinear_gaussian_data_covariance specified.")
        }

        #blocks[["nonlinear_gaussian_data_covariance"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("nonlinear_gaussian_data_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "nonlinear_gaussian_data_variable")
      for (j in 1:length(which_index))
      {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
        }
      }
    }

    if ("data_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "data_variable")
      for (j in 1:length(which_index))
      {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
        }
      }
    }

  }

  for (i in 1:length(blocks[["importance_proposal"]]))
  {
    current_factor = blocks[["importance_proposal"]][[i]]

    if ("evaluate_log_importance_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_importance_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid importance proposal specified.")
        }

        #blocks[["importance_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("simulate_importance_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_importance_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
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
  }

  for (i in 1:length(blocks[["mh_proposal"]]))
  {
    current_factor = blocks[["mh_proposal"]][[i]]

    if ("evaluate_log_mh_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_mh_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid mh proposal specified.")
        }

        blocks[["mh_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("simulate_mh_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_mh_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&","const Data&")) }
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

  }

  for (i in 1:length(blocks[["unadjusted_proposal"]]))
  {
    current_factor = blocks[["unadjusted_proposal"]][[i]]

    if ("simulate_unadjusted_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_unadjusted_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid unadjusted proposal specified.")
        }

      }
      }
    }

  }

  for (i in 1:length(blocks[["independent_mh_proposal"]]))
  {
    current_factor = blocks[["independent_mh_proposal"]][[i]]

    if ("evaluate_log_independent_mh_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_independent_mh_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid independent_mh proposal specified.")
        }

        blocks[["independent_mh_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("simulate_independent_mh_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_independent_mh_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
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

  }

  for (i in 1:length(blocks[["m_proposal"]]))
  {
    current_factor = blocks[["m_proposal"]][[i]]

    if ("simulate_m_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_m_proposal")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE,TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[3] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Parameters&","const Data&")) }
                  , error = function(e) {proposal_type[4] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid m proposal specified.")
        }

        blocks[["m_proposal"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

  }

  for (i in 1:length(blocks[["transition_model"]]))
  {
    current_factor = blocks[["transition_model"]][[i]]

    if ("linear_gaussian_transition_matrix" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_transition_matrix")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c()) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid linear_gaussian_transition_matrix specified.")
        }

        if ("type" %in% names(blocks[["transition_model"]][[i]]))
        {
          if (blocks[["transition_model"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("linear_gaussian_transition_matrix and linear_gaussian_transition_covariance functions are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both linear_gaussian_transition_matrix and linear_gaussian_transition_covariance.")
        }

        #blocks[["linear_gaussian_transition_matrix"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("linear_gaussian_transition_covariance" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_transition_covariance")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        proposal_type = c(TRUE,TRUE)

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c()) }
                  , error = function(e) {proposal_type[1] <<- FALSE})

        tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&")) }
                  , error = function(e) {proposal_type[2] <<- FALSE})

        if (length(which(proposal_type==TRUE))==0)
        {
          stop("No valid linear_gaussian_transition_covariance specified.")
        }

        if ("type" %in% names(blocks[["transition_model"]][[i]]))
        {
          if (blocks[["transition_model"]][[i]][["type"]]!=which(proposal_type)[1])
          {
            stop("linear_gaussian_transition_matrix and linear_gaussian_transition_covariance functions are incompatible.")
          }
        }
        else
        {
          stop("Need to specify both linear_gaussian_transition_matrix and linear_gaussian_transition_covariance.")
        }

        #blocks[["linear_gaussian_transition_covariance"]][[i]][["type"]] = which(proposal_type)[1]
      }
      }
    }

    if ("linear_gaussian_transition_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "linear_gaussian_transition_variable")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
      }
      }
    }

    if ("nonlinear_gaussian_transition_covariance" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "nonlinear_gaussian_transition_covariance")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "arma::mat", c("const Parameters&"))
      }
      }
    }

    if ("nonlinear_gaussian_transition_variable" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "nonlinear_gaussian_transition_variable")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "std::string", c(""))
      }
      }
    }

  }

  for (i in 1:length(blocks[["transition_proposal"]]))
  {
    current_factor = blocks[["transition_proposal"]][[i]]

    if ("evaluate_log_transition_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_transition_proposal")
      for (j in 1:length(which_index))
      {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Data&"))
        }
      }
    }

    if ("simulate_transition_proposal" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "simulate_transition_proposal")
      for (j in 1:length(which_index))
      {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&"))
        }
      }
    }

  }

  for (i in 1:length(blocks[["enk_transform"]]))
  {
    current_factor = blocks[["enk_transform"]][[i]]

    if ("enk_transform" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "enk_transform")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "Parameters", c("const Parameters&"))
      }
      }
    }

    if ("enk_inverse_transform" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "enk_inverse_transform")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "Parameters", c("const Parameters&"))
      }
      }
    }

  }

  for (i in 1:length(blocks[["potential_function"]]))
  {
    current_factor = blocks[["potential_function"]][[i]]

    if ("evaluate_log_potential_function" %in% names(current_factor))
    {
      which_index = which(names(current_factor) == "evaluate_log_potential_function")
      for (j in 1:length(which_index))
      {
      if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
      {
        RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]],  "double", c("const Parameters&","const Parameters&","const Data&"))
      }
      }
    }

  }

  if ("method" %in% names(blocks))
  {
    for (i in 1:length(blocks[["method"]]))
    {
      current_factor = blocks[["method"]][[i]]

      if ("mcmc_weights" %in% names(blocks[["method"]][[i]]))
      {
        which_index = which(names(current_factor) == "mcmc_weights")
        for (j in 1:length(which_index))
        {
        if (inherits(current_factor[[paste("verify_",current_factor[[1]],sep="")]],"XPtr"))
        {
          RcppXPtrUtils::checkXPtr(current_factor[[which_index[j]+1]], "NumericVector")
        }
        }
      }
    }
  }

  return(blocks)
}

strip_out_xptr_stuff <- function(blocks)
{
  if (length(blocks)>0)
  {
    for (i in 1:length(blocks))
    {
      if (length(blocks[[i]])>0)
      {
        for (j in 1:length(blocks[[i]]))
        {
          which_xptr = c()
          current_names = names(blocks[[i]][[j]])
          for (k in 1:length(current_names))
          {
            if (!is.null(current_names))
            {
              if ( grepl("XPtr", current_names[k]) )
              {
                which_xptr = c(which_xptr,k)
              }
            }
          }
          if (length(which_xptr)>0)
          {
            blocks[[i]][[j]] = blocks[[i]][[j]][-which_xptr]
          }
        }
      }
    }
  }

  return(blocks)
}

#' Parse .ilike file to give ilike model.
#'
#' @param filename The name (and path) of the .ilike file containing the model.
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param R_functions (optional) If TRUE, returns R functions. If FALSE, returns C++ functions.
#' @param external_packages (optional) A vector of names of other R packages the functions rely on.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia package that will be installed and loaded.
#' @param verify_cpp_function_types (optional) If TRUE, check the types of the parameters of user-defined .cpp functions. If FALSE (default), types are not checked.
#' @param keep_temporary_model_code (optional) If FALSE (default), the .cpp file generated for compilation is deleted. If TRUE, this file is left in the working directory.
#' @param nesting_level (optional) The level of nesting of the current call to compile. A user should always use the default of 1.
#' @return A list containing the model details.
#' @export
compile <- function(filename,
                    model_parameter_list = list(),
                    R_functions = FALSE,
                    external_packages = c(),
                    julia_bin_dir="",
                    julia_required_libraries=c(),
                    verify_cpp_function_types=FALSE,
                    keep_temporary_model_code=FALSE,
                    nesting_level=1)
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

  model_for_compilation_name = paste("model_for_compilation",nesting_level,"_",ceiling(stats::runif(1, 0, 10^7)),".cpp",sep="")
  fileConn<-file(model_for_compilation_name,open="w")

  the_file = file(filename,open="r")

  in_block = FALSE
  blocks = list(order_of_mcmc=c())
  block_code = ""
  starting_block_flag = FALSE
  pre_blocks = TRUE
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
  unadjusted_proposal_number = 0
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
      blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,model_parameter_list,R_functions,external_packages,julia_bin_dir,julia_required_libraries,fileConn,verify_cpp_function_types,nesting_level,keep_temporary_model_code)
      in_block = FALSE
      starting_block_flag = TRUE
      is_custom = FALSE
      block_code = ""
    }
    else if (nchar(line)>=8)
    {
      if (substr(line, 1, 4)=="/***")
      {

        if (substr(line, nchar(line) - 4 + 1, nchar(line))=="***/")
        {
          # Write initial lines to file.
          if (pre_blocks==TRUE)
          {
            writeLines(c('#include <RcppArmadillo.h>',
              '// [[Rcpp::depends(RcppArmadillo)]]',
              '// [[Rcpp::depends(ilike)]]',
              '// [[Rcpp::depends(BH)]]',
              '// [[Rcpp::depends(dqrng)]]',
              '// [[Rcpp::depends(sitmo)]]',
              '\n',
              'using namespace Rcpp;',
              '#include <ilike.h>',
              'using namespace ilike;',
              '\n'),
              fileConn)

            writeLines(block_code,fileConn)

            pre_blocks = FALSE
            block_code = ""
          }

          # legitimate new section
          starting_block_flag = TRUE

          if (in_block==FALSE) # first block, or previous block wrapped up using //#
          {
            in_block = TRUE
          }
          else # end current block
          {
            blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,model_parameter_list,R_functions,external_packages,julia_bin_dir,julia_required_libraries,fileConn,verify_cpp_function_types,nesting_level,keep_temporary_model_code)
            is_custom = FALSE
            block_code = ""
          }

          unparsed_block_name = substr(line, 5, nchar(line)-4)

          if (nchar(unparsed_block_name)==8)
          {
            stop(paste("Invalid file: line ",line_counter,", new section of file needs a name: use /***name***/.",sep=""))
          }
          split_block_name = split_string_at_comma_ignoring_parentheses(unparsed_block_name)
          new_block_info = determine_block_type(split_block_name,blocks,line_counter,block_type,block_name,factor_number,transition_model_number,potential_function_number,importance_proposal_number,mh_proposal_number,unadjusted_proposal_number,independent_mh_proposal_number,m_proposal_number,enk_transform_number,transition_proposal_number,data_number,method_number)

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
          unadjusted_proposal_number = new_block_info[[11]]
          independent_mh_proposal_number = new_block_info[[12]]
          m_proposal_number = new_block_info[[13]]
          enk_transform_number = new_block_info[[14]]
          transition_proposal_number = new_block_info[[15]]
          data_number = new_block_info[[16]]
          #reinforce_gradient_number = new_block_info[[17]]
          method_number = new_block_info[[17]]

        }
      }
    }

    if ( (pre_blocks) || ( (in_block==TRUE) && (starting_block_flag==FALSE) && (is_custom==TRUE) ) )
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
      blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,model_parameter_list,R_functions,external_packages,julia_bin_dir,julia_required_libraries,fileConn,verify_cpp_function_types,nesting_level,keep_temporary_model_code)
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
      if ( (unadjusted_proposal_number!=0) && unadjusted_proposal_number==length(blocks[["unadjusted_proposal"]]))
      {
        print_unadjusted_proposal_info(length(blocks[["unadjusted_proposal"]]),blocks,line_counter)
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

  if (verify_cpp_function_types)
  {
    blocks = check_types(blocks)
  }

  close(fileConn)

  if (file.exists(model_for_compilation_name))
  {
    Rcpp::sourceCpp(model_for_compilation_name)
  }

  if (keep_temporary_model_code==FALSE)
  {
    file.remove(model_for_compilation_name)
  }

  for (i in 1:length(blocks))
  {
    for (j in 1:length(blocks[[i]]))
    {
      if (length(blocks[[i]][[j]])>0)
      {
        for (k in 1:length(blocks[[i]][[j]]))
        {
          if ( is.character(blocks[[i]][[j]][[k]]) && grepl("XPtr",blocks[[i]][[j]][[k]]) )
          {
            blocks[[i]][[j]][[k+1]] = get(blocks[[i]][[j]][[k]])()
          }
        }
      }
    }
  }

  blocks = strip_out_xptr_stuff(blocks)

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
