#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param for_ggmcmc (optional) Output in tidy format for plotting in ggmcmc package.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            for_ggmcmc = TRUE)
{
  results_directory = paste(results_directory,"/ilike_smc",sep="")

  # Throw error if directory does not exist.
  if (!dir.exists(results_directory))
  {
    stop(paste("Results directory ",results_directory," does not exist.",sep=""))
  }

  # Make the output names.

  names_file = file(paste(results_directory,"/vector_variables.txt",sep=""),open="r")
  line = ""
  while (TRUE)
  {
    previous_line = line
    line = readLines(names_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(names_file)
  variable_names = strsplit(previous_line,";")[[1]]

  sizes_file = file(paste(results_directory,"/vector_variable_sizes.txt",sep=""),open="r")
  while (TRUE)
  {
    previous_line = line
    line = readLines(sizes_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(sizes_file)
  variable_sizes = as.numeric(strsplit(previous_line,";")[[1]])

  lengths_file = file(paste(results_directory,"/output_lengths.txt",sep=""),open="r")
  while (TRUE)
  {
    previous_line = line
    line = readLines(lengths_file, n = 1)
    if ( length(line) == 0 )
    {
      break
    }
  }
  close(lengths_file)
  output_lengths = as.numeric(strsplit(previous_line,split=" +")[[1]])
  output_lengths = output_lengths[!is.na(output_lengths)]
  number_of_chains = length(output_lengths)

  if (length(variable_names)!=length(variable_sizes))
  {
    stop("Variable names and sizes are of different lengths.")
  }

  output_names = rep("",sum(variable_sizes))
  counter = 1
  for (i in 1:length(variable_sizes))
  {
    for (j in 1:variable_sizes[i])
    {
      output_names[counter] = paste(variable_names[i],j,sep="")
      counter = counter + 1
    }
  }

  # Find the final iteration of the SMC algorithm in which the MCMC is stored - this is the folder we need to look in.
  counter = 0
  terminate = FALSE
  iteration_directory = ""
  while (terminate==FALSE)
  {
    previous_iteration_directory = iteration_directory
    iteration_directory = paste(results_directory,"/iteration",counter,sep="")
    if (!dir.exists(iteration_directory))
    {
      terminate = TRUE
    }
    else
    {
      counter = counter + 1
    }
  }

  # Store the output in a data frame.

  # Get number of lines.
  points_filename = paste(previous_iteration_directory,"/vector_points.txt",sep="")
  # points_file = file(points_filename,open="r")

  # # Might not need this part if we write the dimensions to an additional file (see output_lengths file, for example).
  # number_of_lines = 0
  # while (TRUE)
  # {
  #   line = readLines(sizes_file, n = 1)
  #
  #   if (number_of_lines==0)
  #   {
  #     output_cols = length(strsplit(line,",")[[1]])
  #     break
  #   }
  #
  #   # if ( length(line) == 0 )
  #   # {
  #   #   break
  #   # }
  #   # else
  #   # {
  #   #   number_of_lines = number_of_lines + 1
  #   # }
  # }
  # close(points_file)

  all_output_rows = floor(output_lengths)

  if (max(all_output_rows)>0)
  {
    output = read.table(file=points_filename,header=FALSE,sep=",")

    if (nrow(output)!=sum(output_lengths))
    {
      stop("Number of rows in vector_points.txt file does not correspond to output_lengths.txt file.")
    }

    iterations_column = matrix(0,nrow(output),1)
    chains_column = matrix(0,nrow(output),1)

    index = 1
    for (i in 1:length(output_lengths))
    {
      for (j in 1:output_lengths[i])
      {
        iterations_column[index] = j
        chains_column[index] = i
        index = index + 1
      }
    }

    output = cbind(iterations_column,chains_column,output)

    colnames(output) = c('Iteration','Chain',output_names)

    # points_file = file(points_filename,open="r")
    # line_counter = 0
    # point_index = 0
    # chain_index = 1
    # iteration = 1
    #
    # while (TRUE)
    # {
    #   line = readLines(sizes_file, n = 1)
    #   line_counter = line_counter + 1
    #
    #   if ( length(line) == 0 )
    #   {
    #     break
    #   }
    #   else
    #   {
    #     if ((line_counter %% thinning)==0)
    #     {
    #       point_index = point_index + 1
    #       output[point_index,] = c(iteration,chain_index,as.numeric(strsplit(line,",")[[1]]))
    #       iteration = iteration + 1
    #     }
    #   }
    #
    #   if ((line_counter %% sum(output_lengths[1:chain_index]))==0)
    #   {
    #     chain_index = chain_index + 1
    #     iteration = 1
    #   }
    # }
    # close(points_file)

    if (for_ggmcmc==TRUE)
    {
      nParameters =
        output = tidyr::pivot_longer(output, output_names, names_to = "Parameter", values_to = "value")
      attr(output,"nChains") = length(output_lengths)
      attr(output,"nParameters") = length(output_names)
      attr(output,"nIterations") = max(all_output_rows)
      attr(output,"nBurnin") = 0
      attr(output,"nThin") = 1
      attr(output,"description") = "Test"
    }

    return(output)
  }
  else
  {
    stop("No rows found for output.")
  }
}
