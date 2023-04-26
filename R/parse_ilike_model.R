#' Parse .cpp file to give ilike model.
#'
#' @param filename The name (and path) of the .cpp file containing the model.
#' @return A list containing the model details.
#' @export
parse_like_model <- function(filename)
{
  filename = "/Users/richard/Dropbox/projects/ilikemodels/gaussian_unknown_precision/gaussian_unknown_precision_model.cpp"

  the_file = file(filename,open="r")

  in_block = FALSE
  blocks = list()
  block_code = ""
  starting_block_flag = FALSE
  is_custom = FALSE
  line_counter = 0
  add_to_block_code = TRUE

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
            # ignore block if block number is not positive
            if (block_number>0)
            {
              if (is_custom==TRUE)
              {
                blocks[[block_name]][[block_number]] = cppXPtr(block_code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo"))
              }
              else
              {
                blocks[[block_name]][[block_number]][["type"]] = block_type
                blocks[[block_name]][[block_number]][["variables"]] = variables
                blocks[[block_name]][[block_number]][["parameters"]] = parameters
              }
            }
            block_code = ""
          }

          if (nchar(line)==8)
          {
            stop("New section of file needs a name: use /***name***/.")
          }

          unparsed_block_name = substr(line, 5, nchar(line)-4)

          split_block_name = strsplit(unparsed_block_name,split=",")[[1]]

          if (length(split_block_name)==1)
          {
            block_name = split_block_name
            block_number = 1
            is_custom = TRUE
          }
          else if (length(split_block_name)==2)
          {
            block_number = strtoi(split_block_name[2])
            if (is.na(block_number))
            {
              stop(paste("Invalid file: line ",line_number,", expected number as second argument.",sep=""))
            }
            else
            {
              is_custom = TRUE
            }
            block_name = split_block_name[1]
          }
          else if (length(split_block_name)>=3)
          {
            is_custom = FALSE
            block_name = split_block_name[1]

            block_number = strtoi(split_block_name[2])
            if (is.na(block_number))
            {
              block_number = 1
              block_type_index = 2
              #stop(paste("Invalid file: line ",line_number,", expected number as second argument.",sep=""))
            }
            else
            {
              block_type_index = 3
              if (length(split_block_name)==3)
              {
                stop(paste("Invalid file: line ",line_number,", require variables and parameters to be specified.",sep=""))
              }
            }

            block_type = split_block_name[block_type_index]

            variables = split_block_name[block_type_index+1]

            parameters = list()
            if (block_type_index+1<length(split_block_name))
            {
              for (i in 2:(length(split_block_name)-block_type_index))
              {
                index = block_type_index + i
                parameters[[i-1]] = split_block_name[index]
              }
            }
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
      if (is_custom==TRUE)
      {
        blocks[[block_name]][[block_number]] = RcppXPtrUtils::cppXPtr(block_code,plugins=c("cpp11"),depends = c("ilike","RcppArmadillo","BH","dqrng","sitmo"))
      }
      else
      {
        blocks[[block_name]][[block_number]][["type"]] = block_type
        blocks[[block_name]][[block_number]][["variables"]] = variables
        blocks[[block_name]][[block_number]][["parameters"]] = parameters
      }
    }
  }

  close(the_file)

  if ("evaluate_log_prior" %in% names(blocks))
  {
    for (i in 1:length(blocks$evaluate_log_prior))
    {
      RcppXPtrUtils::checkXPtr(blocks$evaluate_log_prior[[i]], "double", c("const Parameters&"))
    }
  }

  if ("simulate_prior" %in% names(blocks))
  {
    for (i in 1:length(blocks$simulate_prior))
    {
      RcppXPtrUtils::checkXPtr(blocks$simulate_prior[[i]], "Parameters", c("RandomNumberGenerator&"))
    }
  }

  if ("evaluate_log_likelihood" %in% names(blocks))
  {
    for (i in 1:length(blocks$evaluate_log__likelihood))
    {
      RcppXPtrUtils::checkXPtr(blocks$evaluate_log__likelihood[[i]], "double", c("const Parameters&","const Data&"))
    }
  }

  return(blocks)
}

# expect input for each block in one of the following forms:
# (a) /***evaluate_log_prior***/, followed by a C++ function
# (b) /***evaluate_log_prior,n***/, in the case of a model with multiple factors, where this is the nth factor
# (c) /***evaluate_log_prior,norm,0,1***/, to use a normal distribution with parameters 0 and 1
# (d) /***evaluate_log_prior,n,norm,0,1***/, to use a normal distribution with parameters 0 and 1, in the case of a model with multiple factors, where this is the nth factor
