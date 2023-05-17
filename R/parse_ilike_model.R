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

extract_block <- function(blocks,block_name,block_number,block_code,block_function,is_custom,parameter_list)
{
  #browser()

  # ignore block if block number is not positive
  if (block_number>0)
  {
    if (is_custom==TRUE)
    {
      blocks[[block_name]][[block_number]] = RcppXPtrUtils::cppXPtr(block_code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo"))
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
          stop(paste("Block ",block_name,", number ",block_number,": input does not specify a function.",sep=""))
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
        stop(paste("Block ",block_name,", number ",block_number,": function needs to have a name.",sep=""))
      }

      if (is_like_function==TRUE)
      {
        if (block_name=="prior")
        {
          split_arg_string = strsplit(arg_string,",")

          variables = split_arg_string[1]

          parameters = list()
          if (length(split_arg_string)>1)
          for (k in 2:(length(split_arg_string)))
          {
            parameters[[k-1]] = split_arg_string[k]
          }

          blocks[[block_name]][[block_number]][["type"]] = ilike_type
          blocks[[block_name]][[block_number]][["variables"]] = variables
          blocks[[block_name]][[block_number]][["parameters"]] = parameters
        }
        else
        {
          stop(paste("Block ",block_name,", number ",block_number,": ",block_name," is invalid block type.",sep=""))
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
          stop(paste("Block ",block_name,", number ",block_number,": can have a maximum of one equals sign in it.",sep=""))
        }

        function_info = ilike_parse(function_to_use,
                                    parameter_list)

        R_function_arguments = function_info[[2]]

        cpp_function_arguments_string = ""
        which_parameters = c()
        which_data = c()
        if (length(R_function_arguments)>0)
        {
          which_parameters = matrix(0,length(R_function_arguments))
          which_data = matrix(0,length(R_function_arguments))

          for (i in 1:length(R_function_arguments))
          {
            split_at_dot = strsplit(R_function_arguments[i],"\\.")[[1]]
            if (length(split_at_dot)==0)
            {
              stop(paste("Block ",block_name,", number ",block_number,": argument is of size zero.",sep=""))
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

        cpp_function_name = paste(block_name,block_number,sep="")
        R_function_name = paste(cpp_function_name,'_R',sep="")

        temp_filename = paste(cpp_function_name,".cpp",sep="")
        xptr_name = paste(cpp_function_name,"getXPtr",sep="_")

        if (block_name=="data")
        {
          if (sum(which_parameters)>0)
          {
            stop(paste("Block ",block_name,", number ",block_number,": using variables from parameters not possible in a data function.",sep=""))
          }

          return_type = "Data"
          arguments = c()
          #args_for_typedef = ""
          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", number ",block_number,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Data output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else if (block_name=="evaluate_log_prior")
        {
          if (sum(which_data)>0)
          {
            stop(paste("Block ",block_name,", number ",block_number,": using variables from data not possible in a prior.",sep=""))
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
            stop(paste("Block ",block_name,", number ",block_number,": using variables from data not possible in a prior.",sep=""))
          }

          if (output_variable=="")
          {
            stop(paste("Block ",block_name,", number ",block_number,": no output variables specified, need output_variable=some_function(...).",sep=""))
          }

          return_type = "Parameters"
          arguments = c()
          arguments[1] = "RandomNumberGenerator &rng"
          #args_for_typedef = "RandomNumberGenerator&"
          function_body = paste('Function f(Environment::global_env()["',R_function_name,'"]); Parameters output; output["',output_variable,'"] = NumericVector(f(',cpp_function_arguments_string,')); return output;',sep="")
        }
        else
        {
          stop(paste("Block ",block_name,", number ",block_number,": block is of unknown type.",sep=""))
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

        blocks[[block_name]][[block_number]] = get(xptr_name)()
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
            blocks = extract_block(blocks,block_name,block_number,block_code,block_function,is_custom,parameter_list)
            is_custom = FALSE
            block_code = ""
          }

          if (nchar(line)==8)
          {
            stop("New section of file needs a name: use /***name***/.")
          }

          unparsed_block_name = substr(line, 5, nchar(line)-4)

          # expect input for each block in one of the following forms:
          # (a) /***evaluate_log_prior***/, followed by a C++ function
          # (b) /***evaluate_log_prior,n***/, in the case of a model with multiple factors, where this is the nth factor
          # (c) /***prior,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1
          # (d) /***prior,n,ilike::lnorm(tau,1,1)***/, to use a lognormal distribution with parameters 1 and 1, in the case of a model with multiple factors, where this is the nth factor
          # (e) /***evaluate_log_prior,dlnorm(tau,1,1,TRUE)***/, so that we extract the named argument from the parameters read into the function, and input the other arguments, in this order, in to dlnorm (dlnorm in this case is a base R function, but it could be a user defined one in an R, python or julia file)
          # (f) /***evaluate_log_prior,n,dlnorm(tau,1,1,TRUE)***/, is the same case as (e), except that we have a model with multiple factors and this is the nth factor
          # (g) /***evaluate_log_likelihood,n,dnorm(data.y,0,recip(sqrt(tau)),TRUE)***/
          # (h) /***simulate_prior,tau=rlnorm(1,1,1)***/, when we simulate a variable
          # (i) /***simulate_prior,1,theta=rlnorm(1,theta,1)***/, for simulations that involve multiple steps
          # (j) /***simulate_prior,2,tau=rlnorm(1,theta,1)***/, the counterpart to (i)

          split_block_name = split_string(unparsed_block_name)

          if (length(split_block_name)==1)
          {
            block_name = split_block_name
            block_number = 1
            is_custom = TRUE
          }
          else if (length(split_block_name)==2)
          {
            # block number is always second argument, if present
            block_number = strtoi(split_block_name[2])
            if (is.na(block_number)) # not a number, so second argument must be a function
            {
              block_function = split_block_name[2]
              block_number = 1
            }
            else # is a number. There are no other arguments, so we have a custom C++ function.
            {
              is_custom = TRUE
            }
            block_name = split_block_name[1]
          }
          else if (length(split_block_name)==3)
          {
            is_custom = FALSE
            block_name = split_block_name[1]

            block_number = strtoi(split_block_name[2])
            if (is.na(block_number))
            {
              # only way of having 3 arguments is that the second is a number
              stop(paste("Invalid file: line ",line_counter,", expected number as second argument.",sep=""))
            }
            else
            {
              block_function = split_block_name[3]
            }
          }
          else
          {
            stop(paste("Invalid file: line ",line_counter,", invalid function definition (can only have up three parts, separated by commas).",sep=""))
          }
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
    if (block_number>0)
    {
      blocks = extract_block(blocks,block_name,block_number,block_code,block_function,is_custom,parameter_list)
    }
  }

  close(the_file)

  # if ("data" %in% names(blocks))
  # {
  #   for (i in 1:length(blocks$data))
  #   {
  #     RcppXPtrUtils::checkXPtr(blocks$data[[i]], "Data")
  #   }
  # }
  #
  # if ("evaluate_log_prior" %in% names(blocks))
  # {
  #   for (i in 1:length(blocks$evaluate_log_prior))
  #   {
  #     RcppXPtrUtils::checkXPtr(blocks$evaluate_log_prior[[i]], "double", c("const Parameters&"))
  #   }
  # }
  #
  # if ("simulate_prior" %in% names(blocks))
  # {
  #   for (i in 1:length(blocks$simulate_prior))
  #   {
  #     RcppXPtrUtils::checkXPtr(blocks$simulate_prior[[i]], "Parameters", c("RandomNumberGenerator&"))
  #   }
  # }
  #
  # if ("evaluate_log_likelihood" %in% names(blocks))
  # {
  #   for (i in 1:length(blocks$evaluate_log_likelihood))
  #   {
  #     RcppXPtrUtils::checkXPtr(blocks$evaluate_log_likelihood[[i]], "double", c("const Parameters&","const Data&"))
  #   }
  # }

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
