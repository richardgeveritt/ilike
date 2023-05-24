split_string <- function(input_string) {
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
          stop(paste("In call ",arg_string,", parameter ",as.numeric(substr(h,2,nchar(h)))," not found in parameter_list (parameter_list has length ",length(parameter_list),".",sep=""))
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

factor_processing = function(factor_number,blocks,block_name,prior_function_types,line_counter)
{
  # Is this a continuation of the current factor, or a new one?

  # Get the current factor info.
  if ("factor" %in% names(blocks))
  {
    current_factor_info = blocks[["factor"]][[factor_number]]

    # Possible factors are:
    # (a) "prior" (prior variant)
    # (b) one or both of "evaluate_log_prior","simulate_prior" (prior variant)
    # (c) "evaluate_log_likelihood", alone (likelihood variant)
    # (d) "simulate_model" and (ABC, SL, etc) (likelihood variant)

    if (block_name %in% prior_function_types)
    {
      new_is_llhd = FALSE
    }
    else
    {
      new_is_llhd = TRUE
    }

    if ("prior" %in% names(current_factor_info))
    {
      # factor is complete
      print_factor_info(factor_number,blocks,line_counter-1)
      factor_number = factor_number + 1
    }
    else if ( ("evaluate_log_prior" %in% names(current_factor_info)) || ("simulate_prior" %in% names(current_factor_info)) )
    {
      # factor is complete
      if ( ("evaluate_log_prior" %in% names(current_factor_info)) && (block_name=="evaluate_log_prior") )
      {
        # if ("simulate_prior" %in% names(current_factor_info))
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_prior and simulate_prior.',sep = ""))
        # }
        # else
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_prior.',sep = ""))
        # }
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if ( ("simulate_prior" %in% names(current_factor_info)) && (block_name=="simulate_prior") )
      {
        # if ("evaluate_log_prior" %in% names(current_factor_info))
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_prior and simulate_prior.',sep = ""))
        # }
        # else
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains simulate_prior.',sep = ""))
        # }
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if (new_is_llhd)
      {
        # if ( ("evaluate_log_prior" %in% names(current_factor_info)) && ("simulate_prior" %in% names(current_factor_info)) )
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_prior and simulate_prior.',sep = ""))
        # }
        # else if ("evaluate_log_prior" %in% names(current_factor_info))
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_prior.',sep = ""))
        # }
        # else if ("simulate_prior" %in% names(current_factor_info))
        # {
        #   print(paste('Factor ends on line ',line_counter-1,'. Contains simulate_prior.',sep = ""))
        # }
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
      else if (block_name=="prior")
      {
        print_factor_info(factor_number,blocks,line_counter-1)
        factor_number = factor_number + 1
      }
    }
    else if ("evaluate_log_likelihood" %in% names(current_factor_info))
    {
      # factor is complete
      #print(paste('Factor ends on line ',line_counter-1,'. Contains evaluate_log_likelihood.',sep = ""))
      print_factor_info(factor_number,blocks,line_counter-1)
      factor_number = factor_number + 1
    }
  }
  else
  {
    factor_number = factor_number + 1
  }

  return(factor_number)
}

data_processing = function(data_number)
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
        stop(paste("Invalid file: line ",line_counter,', importance_proposal specified, but previous importance proposal block is incomplete (did you specify both evaluate_log_importance_proposal and simulate_importance_proposal?',sep=""))
      }

    }
  }
  else
  {
    importance_proposal_number = importance_proposal_number + 1
  }

  return(importance_proposal_number)
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

determine_block_type = function(split_block_name,blocks,line_counter,block_type,block_name,factor_number,importance_proposal_number,data_number)
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

  prior_function_types = c("prior","evaluate_log_prior","simulate_prior")
  likelihood_function_types = c("evaluate_log_likelihood","simulate_model")
  factor_function_types = c(prior_function_types,likelihood_function_types)
  data_function_types = c("data")
  importance_proposal_types = c("simulate_importance_proposal","evaluate_log_importance_proposal")

  # distinguish between factor, data, etc
  if (block_name %in% factor_function_types)
  {
    block_type = "factor"
    factor_number = factor_processing(factor_number,blocks,block_name,prior_function_types,line_counter)
    number_to_pass_to_extract_block = factor_number
  }
  else if (block_name %in% data_function_types)
  {
    block_type = "data"
    data_number = data_processing(data_number)
    number_to_pass_to_extract_block = data_number
  }
  else if (block_name %in% importance_proposal_types)
  {
    block_type = "importance_proposal"
    importance_proposal_number = importance_proposal_processing(importance_proposal_number,blocks,block_name,line_counter)
    number_to_pass_to_extract_block = importance_proposal_number
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
              importance_proposal_number,
              data_number))
}

extract_block <- function(blocks,block_type,block_name,factor_number,line_counter,block_code,block_function,is_custom,parameter_list)
{
  # ignore block if block number is not positive
  if (factor_number>0)
  {
    if (is_custom==TRUE)
    {
      my_list = list(RcppXPtrUtils::cppXPtr(block_code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo")))
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
        if (block_name=="prior")
        {
          split_arg_string = strsplit(arg_string,",")[[1]]

          variables = split_arg_string[1]

          parameters = list()
          if (length(split_arg_string)>1)
          for (k in 2:(length(split_arg_string)))
          {
            parameters[[k-1]] = split_arg_string[k]
          }

          # temp_list = blocks[[block_type]]
          # temp_list[[factor_number]][[block_name]][["type"]] = ilike_type
          # temp_list[[factor_number]][[block_name]][["variables"]] = variables
          # temp_list[[factor_number]][[block_name]][["parameters"]] = parameters
          # blocks[[block_type]] = temp_list

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

          #blocks[[block_type]][[factor_number]][[block_name]][["type"]] = ilike_type
          #blocks[[block_type]][[factor_number]][[block_name]][["variables"]] = variables
          #blocks[[block_type]][[factor_number]][[block_name]][["parameters"]] = parameters
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
        which_parameters = c()
        which_proposal_parameters = c()
        which_data = c()
        if (length(R_function_arguments)>0)
        {
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

        if (block_name=="data")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": using variables from parameters not possible in a data function.",sep=""))
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

          return_type = "double"
          arguments = c()
          arguments[1] = "const Parameters &parameters"
          #args_for_typedef = "const Parameters&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="evaluate_log_likelihood")
        {
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
          if ( (sum(which_data)>0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            arguments[3] = "const Data &data"
            #args_for_typedef = "const Parameters&"
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"
          }
          else
          {
            return_type = "double"
            arguments = c()
            arguments[1] = "const Parameters &parameters"
            #args_for_typedef = "const Parameters&"
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); ','return NumericVector(f(',cpp_function_arguments_string,'))[0];',sep="")
        }
        else if (block_name=="simulate_importance_proposal")
        {
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", line number ",line_counter,": no output variables specified, need output_variable=some_function(...).",sep=""))
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
          }
          else if ( (sum(which_data)>0) && (sum(which_proposal_parameters)==0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Data &data"
            #args_for_typedef = "const Parameters&"
          }
          else if ( (sum(which_data)==0) && (sum(which_proposal_parameters)>0) )
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            arguments[2] = "const Parameters &proposal_parameters"
            #args_for_typedef = "const Parameters&"
          }
          else
          {
            arguments = c()
            arguments[1] = "RandomNumberGenerator &rng"
            #args_for_typedef = "const Parameters&"
          }

          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
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
        close(fileConn)

        Rcpp::sourceCpp(temp_filename)
        file.remove(temp_filename)

        my_list = list()
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
          blocks[[block_type]][[factor_number]] = append(blocks[[block_type]] [[factor_number]],my_list)
        }
        blocks[[block_type]][[factor_number]][[block_name]] = get(xptr_name)
      }

    }
  }

  return(blocks)
}

#' Parse .cpp file to give ilike model.
#'
#' @param filename The name (and path) of the .cpp file containing the model.
#' @param parameter_list (optional) A list containing parameters for the model.
#' @param julia_bin_dir (optional) The directory containing the Julia bin file - only needed if Julia functions are used.
#' @param julia_required_libraries (optional) Vector of strings, each of which is a Julia packge that will be installed and loaded.
#' @return A list containing the model details.
#' @export
parse_ilike_model <- function(filename,
                              parameter_list = list(),
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
  blocks = list()
  block_code = ""
  starting_block_flag = FALSE
  is_custom = FALSE
  line_counter = 0
  add_to_block_code = TRUE
  function_names = c()
  in_factor = FALSE

  factor_number = 0
  importance_proposal_number = 0
  data_number = 0
  block_type = "none"
  block_name = "none"

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

    if (nchar(line)>=8)
    {
      if (substr(line, 1, 4)=="/***")
      {
        if (substr(line, nchar(line) - 4 + 1, nchar(line))=="***/")
        {
          # legitimate new section
          starting_block_flag = TRUE

          if (in_block==FALSE) # first block
          {
            in_block = TRUE
          }
          else # end current block
          {
            blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,parameter_list)
            is_custom = FALSE
            block_code = ""
          }

          unparsed_block_name = substr(line, 5, nchar(line)-4)

          if (nchar(unparsed_block_name)==8)
          {
            stop(paste("Invalid file: line ",line_counter,", new section of file needs a name: use /***name***/.",sep=""))
          }
          split_block_name = split_string(unparsed_block_name)
          new_block_info = determine_block_type(split_block_name,blocks,line_counter,block_type,block_name,factor_number,importance_proposal_number,data_number)

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
          importance_proposal_number = new_block_info[[7]]
          data_number = new_block_info[[8]]

        }
      }
    }

    if ( (in_block==TRUE) && (starting_block_flag==FALSE) && (is_custom==TRUE) )
    {
      block_code = paste(block_code,line,sep="\n")
    }

  }

  if (in_block==TRUE)
  {
    # ignore block if block number is not positive
    if (number_to_pass_to_extract_block>0)
    {
      blocks = extract_block(blocks,block_type,block_name,number_to_pass_to_extract_block,line_counter,block_code,block_function,is_custom,parameter_list)
      if (factor_number==length(blocks[["factor"]]))
      {
        print_factor_info(length(blocks[["factor"]]),blocks,line_counter)
      }
      if (importance_proposal_number==length(blocks[["importance_proposal"]]))
      {
        print_importance_proposal_info(length(blocks[["importance_proposal"]]),blocks,line_counter)
      }
    }
  }

  close(the_file)

  if ("data" %in% names(blocks))
  {
    for (i in 1:length(blocks[["data"]]))
    {
      RcppXPtrUtils::checkXPtr(blocks[["data"]][[i]][["data"]], "Data")
    }
  }

  for (i in 1:length(blocks[["factor"]]))
  {
    current_factor = blocks[["factor"]][[i]]
    if ("evaluate_log_prior" %in% names(current_factor))
    {
      RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_prior"]], "double", c("const Parameters&"))
    }

    if ("simulate_prior" %in% names(current_factor))
    {
      RcppXPtrUtils::checkXPtr(current_factor[["simulate_prior"]], "Parameters", c("RandomNumberGenerator&"))
    }

    if ("evaluate_log_prior" %in% names(current_factor))
    {
      RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_prior"]], "double", c("const Parameters&"))
    }

    if ("evaluate_log_likelihood" %in% names(current_factor))
    {
      RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_likelihood"]],  "double", c("const Parameters&","const Data&"))
    }

  }

  for (i in 1:length(blocks[["importance_proposal"]]))
  {
    current_importance_proposal = blocks[["importance_proposal"]][[i]]

    if ("evaluate_log_importance_proposal" %in% names(current_importance_proposal))
    {
      proposal_type = c(TRUE,TRUE)#,TRUE,TRUE)

      tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&")) }
                , error = function(e) {proposal_type[1] <<- FALSE})

      # tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Parameters&")) }
      #           , error = function(e) {proposal_type[2] <<- FALSE})

      tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Data&")) }
                , error = function(e) {proposal_type[2] <<- FALSE})

      # tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[["evaluate_log_importance_proposal"]],  "double", c("const Parameters&","const Parameters&","const Data&")) }
      #           , error = function(e) {proposal_type[4] <<- FALSE})

      if (length(which(proposal_type==TRUE))==0)
      {
        stop("No valid importance proposal specified.")
      }

      blocks[["importance_proposal"]][[i]][["type"]] = which(proposal_type)[1]
    }

    if ("simulate_importance_proposal" %in% names(current_importance_proposal))
    {
      proposal_type = c(TRUE,TRUE)#TRUE,TRUE)

      tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&")) }
                , error = function(e) {proposal_type[1] <<- FALSE})

      # tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&")) }
      #           , error = function(e) {proposal_type[2] <<- FALSE})

      tryCatch( { RcppXPtrUtils::checkXPtr(current_importance_proposal[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Data&")) }
                , error = function(e) {proposal_type[2] <<- FALSE})

      # tryCatch( { RcppXPtrUtils::checkXPtr(current_factor[["simulate_importance_proposal"]], "Parameters", c("RandomNumberGenerator&","const Parameters&","const Data&")) }
      #           , error = function(e) {proposal_type[4] <<- FALSE})

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
