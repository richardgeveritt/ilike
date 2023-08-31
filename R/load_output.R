#' Loading SMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggsmc (optional) Output in tidy format for plotting in gggsmc package.
#' @param as.mcmc (optional) Output treats particles as different MCMC chains.
#' @param as.enk (optional) Output treats particles as an ensemble.
#' @param which.targets (optional) The indices of the targets to output (defaults to all).
#' @return A list containing the SMC output.
#' @export
load_smc_output = function(results_directory,
                           ggsmc = TRUE,
                           as.mcmc = FALSE,
                           as.enk = FALSE,
                           which.targets = NULL)
{
  if (as.mcmc && as.enk)
  {
    stop("Cannot treat output as both MCMC chains and an ensemble.")
  }

  description = results_directory
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
  output_lengths = vector("list", length(variable_sizes))
  counter = 1
  line = ""
  while (TRUE)
  {
    previous_line = line
    line = readLines(lengths_file, n = 1)

    if ( length(line) == 0 )
    {
      break
    }
    else
    {
      output_lengths[[counter]] = as.numeric(strsplit(line,split=" +")[[1]])
      output_lengths[[counter]] = output_lengths[[counter]][!is.na(output_lengths[[counter]])]

      if (sum(output_lengths[[counter]])==0)
      {
        stop("Output reported to be of length 0 in output_lengths.txt. Error in output files.")
      }

      counter = counter + 1
    }
  }
  close(lengths_file)

  #number_of_points = length(output_lengths)

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

  # Store the output in a data frame.

  # Find the final iteration of the SMC algorithm in which the MCMC is stored - this is the folder we need to look in.
  all_dirs = list.dirs(results_directory,recursive = FALSE)

  if (is.null(which.targets))
  {
    if (!as.mcmc)
    {
      which.targets = 1:length(all_dirs)
    }
    else
    {
      which.targets = length(all_dirs)
    }
  }
  else
  {
    inputted_which_targets = which.targets

    target_indices = matrix(0,length(all_dirs))
    for (k in 1:length(all_dirs))
    {
      split = strsplit(all_dirs[k],"/")[[1]]
      target_indices[i] = strtoi(gsub("iteration","",split[length(split)]))
    }

    which.targets = matrix(0,length(inputted_which_targets))
    for (k in 1:length(inputted_which_targets))
    {
      which.targets = which(target_indices==inputted_which_targets[k])
    }
  }

  for (k in which.targets)
  {
    split = strsplit(all_dirs[k],"/")[[1]]
    target = strtoi(gsub("iteration","",split[length(split)]))

    iteration_directory = all_dirs[k]

    points_filename = paste(iteration_directory,"/vector_points.txt",sep="")

    output = read.table(file=points_filename,header=FALSE,sep=",")

    if (nrow(output)!=sum(output_lengths[[k]]))
    {
      stop("Number of rows in vector_points.txt file does not correspond to output_lengths.txt file.")
    }

    schedule_parameters_filename = paste(iteration_directory,"/schedule_parameters.txt",sep="")
    tryCatch( {schedule_parameters = read.table(file=schedule_parameters_filename,header=FALSE) }
              , error = function(e) {schedule_parameters <<- NULL})
    if (is.null(schedule_parameters))
    {
      schedule_parameters_column = matrix(0,nrow(output),1)
    }
    else
    {
      schedule_parameters_column = matrix(paste(schedule_parameters),nrow(output),1)
    }

    iterations_column = matrix(0,nrow(output),1)
    chains_column = matrix(0,nrow(output),1)
    target_column = matrix(0,nrow(output),1)

    index = 1
    for (i in 1:length(output_lengths[[k]]))
    {
      for (j in 1:output_lengths[[k]][i])
      {
        iterations_column[index] = j
        chains_column[index] = i
        target_column[index] = target
        index = index + 1
      }
    }

    if (as.mcmc)
    {
      output = cbind(iterations_column,chains_column,output)
      colnames(output) = c('Iteration','Chain',output_names)
    }
    else
    {
      if (as.enk)
      {
        if (is.null(schedule_parameters))
        {
          output = cbind(target_column,iterations_column,chains_column,output)
          colnames(output) = c('Target','Iteration','Particle',output_names)
        }
        else
        {
          output = cbind(target_column,schedule_parameters_column,iterations_column,chains_column,output)
          colnames(output) = c('Target','TargetParameters','Iteration','Particle',output_names)
        }
      }
      else
      {
        log_weight_filename = paste(iteration_directory,"/normalised_log_weights.txt",sep="")
        log_weight = read.table(file=log_weight_filename,header=FALSE,sep=",")

        ancestor_index_filename = paste(iteration_directory,"/ancestor_index.txt",sep="")
        tryCatch( {ancestor_index = read.table(file=ancestor_index_filename,header=FALSE,sep=",") + 1 }
                  , error = function(e) {ancestor_index <<- 1:(nrow(log_weight))})

        if (is.null(schedule_parameters))
        {
          output = cbind(target_column,iterations_column,chains_column,ancestor_index,log_weight,output)
          colnames(output) = c('Target','Iteration','Particle','AncestorIndex','LogWeight',output_names)
        }
        else
        {
          output = cbind(target_column,schedule_parameters_column,iterations_column,chains_column,ancestor_index,log_weight,output)
          colnames(output) = c('Target','TargetParameters','Iteration','Particle','AncestorIndex','LogWeight',output_names)
        }
      }
    }

    if ((k==1) || (as.mcmc))
    {
      all_output = output
    }
    else
    {
      all_output = rbind(all_output,output)
    }

  }

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

  if ( (ggsmc==TRUE) && (as.mcmc==TRUE) )
  {
    nParameters =
      all_output = tidyr::pivot_longer(all_output, all_of(output_names), names_to = "Parameter", values_to = "Value")
    attr(output,"nChains") = length(output_lengths[[k]])
    attr(output,"nParameters") = length(output_names)
    attr(output,"nIterations") = max(output_lengths[[k]])
    attr(output,"nBurnin") = 0
    attr(output,"nThin") = 1
    attr(output,"description") = description
  }

  if ( (ggsmc==TRUE) && (as.mcmc==FALSE) )
  {
    nParameters =
      all_output = tidyr::pivot_longer(all_output, all_of(output_names), names_to = "Parameter", values_to = "Value")
    attr(all_output,"nTargets") = length(all_dirs)
    attr(all_output,"nParticles") = length(output_lengths[[k]])
    attr(all_output,"nParameters") = length(output_names)
    attr(all_output,"nIterations") = max(output_lengths[[k]])
    attr(all_output,"nBurnin") = 0
    attr(all_output,"nThin") = 1
    attr(all_output,"description") = description
  }

  all_output$Target = as.integer(all_output$Target)

  return(all_output)
}

#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggmcmc (optional) Output in tidy format for plotting in ggmcmc package.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            ggmcmc = TRUE)
{
  return(load_smc_output(results_directory = results_directory,ggsmc = ggmcmc,as.mcmc=TRUE))
}

#' Loading ensemble Kalman output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggsmc (optional) Output in tidy format for plotting in ggsmc package.
#' @return A list containing the ensemble members (called particles).
#' @export
load_enk_output = function(results_directory,
                            ggsmc = TRUE)
{
  return(load_smc_output(results_directory = results_directory,ggsmc = ggsmc,as.enk=TRUE))
}
