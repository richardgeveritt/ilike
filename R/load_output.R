extract_numbers <- function(inputString)
{
  # Use regular expression to find all numbers in the string
  matches <- gregexpr("\\d+", inputString, perl = TRUE)

  # Extract the matched numbers
  numbers <- regmatches(inputString, matches)[[1]]

  # Convert the character vector to numeric
  numbers <- as.numeric(numbers)

  return(numbers)
}

order_directories <- function(all_dirs)
{
  numbers = matrix(0,length(all_dirs))
  for (k in 1:length(all_dirs))
  {
    split = strsplit(all_dirs[k],"/")[[1]]
    numbers[k] = extract_numbers(split[length(split)])
  }

  return(all_dirs[sort(numbers,index.return=TRUE)$ix])
}

get_directory_numbers <- function(all_dirs)
{
  numbers = matrix(0,length(all_dirs))

  for (k in 1:length(all_dirs))
  {
    split = strsplit(all_dirs[k],"/")[[1]]
    numbers[k] = extract_numbers(split[length(split)])
  }

  return(numbers)
}

#' Take output from the functions load_mcmc_output, load_smc_output or load_enk_output and convert it to a more standard format (non-tidy format with one column per dimension).
#'
#' @param output Output from the functions load_mcmc_output, load_smc_output or load_enk_output.
#' @param variables (optional) Variables to include in the output (default is all).
#' @export
ilike_pivot_wider <- function(output,
                              variables=NULL)
{
  if (!is.null(variables))
  {
    output = poorman::filter(output,ParameterName %in% variables)
  }

  new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },output$ParameterName,output$Dimension)
  output = subset(output,select = -c(ParameterName,Dimension))
  output$Parameter = new_variable_names
  output = poorman::distinct(output)
  return(poorman::pivot_wider(output, names_from = "Parameter", values_from = "Value"))
}

#' Loading SMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggmcmc (optional) Output in tidy format for plotting in ggmcmc package.
#' @param ggsmc (optional) Output in tidy format for plotting in ggsmc package.
#' @param ilike.output (optional) Output can be processed by ilike,output package (default is TRUE). Choosing TRUE for this argument is incompatiable with ggmcmc=TRUE, since the two packages require different formatting of the output.
#' @param as.mcmc (optional) Output treats particles as different MCMC chains.
#' @param as.enk (optional) Output treats particles as an ensemble.
#' @param which.targets (optional) The indices of the targets to output (defaults to all).
#' @param directory_prefix (optional; for nested output only) The first part of the name of the directory within results_directory that contains the results. (default is "ilike", giving a directory of results_directory/ilike_smc)
#' @param external_log_weights (optional; for nested output only) The weights of the importance points external to the current folder. (default is 1, to be used at the top level of nested output)
#' @param external_target_parameters (optional; for nested output only) The parameters of the target external to the current folder. (default is "", corresponding to no parameters)
#' @param nesting_level (optional; for nested output only) The level of nesting at which to extract points. (default is 0, representing the top level of nested output)
#' @param factor (optional; for nested output only) The factor from which to extract points. (default is 1)
#' @return A list containing the SMC output.
#' @export
load_smc_output = function(results_directory,
                           ggmcmc = FALSE,
                           ggsmc = TRUE,
                           ilike.output = TRUE,
                           as.mcmc = FALSE,
                           as.enk = FALSE,
                           which.targets = NULL,
                           directory_prefix = "ilike",
                           external_log_weights = c(0),
                           external_target_parameters = "",
                           nesting_level = 0,
                           factor = 1)
{
  description = results_directory

  if (nesting_level==0)
  {
    if (as.enk==TRUE)
    {
      results_directory = paste(results_directory,"/",directory_prefix,"_enk",sep="")
    }
    else
    {
      results_directory = paste(results_directory,"/",directory_prefix,"_smc",sep="")
    }
  }
  else
  {
    # Can't yet cope with enk nested.
    results_directory = paste(results_directory,"/",directory_prefix,"_smc",sep="")
  }

  # Throw error if directory does not exist.
  if (!dir.exists(results_directory))
  {
    stop(paste("Results directory ",results_directory," does not exist.",sep=""))
  }

  all_dirs = list.dirs(results_directory,recursive = FALSE)
  all_dirs = order_directories(all_dirs)

  # Factor name.
  factor_name = paste("factor",factor,sep="")

  # If nesting_level == 0, then we are at the top level of the nesting hierarchy.
  # We need to run the main body of the function.
  # If nesting_level != 0, we are not at the top level.
  # We need to run load_smc_output at the top level to get the external points.
  # We need to find a list of subdirectories, and call load_smc_output (using external points input)
  # in each subdirectory, and concatenate the output.
  # For now, throw error if nesting_level>1, since then we need to do something additional to make the correct external_weights
  # to pass into the function.
  if (nesting_level==1)
  {
    for (i in 1:length(all_dirs))
    {
      sub_output = load_smc_output(results_directory=all_dirs[i],
                                   ggmcmc=ggmcmc,
                                   ggsmc=ggsmc,
                                   as.mcmc = as.mcmc,
                                   as.enk = as.enk,
                                   which.targets = which.targets,
                                   directory_prefix=factor_name,
                                   external_log_weights = external_log_weights,#poorman::filter(poorman::filter(theta_output,Target==i),Dimension==1)$LogWeight,
                                   nesting_level = nesting_level-1)


      numbers = extract_numbers(all_dirs[i])
      target_id = numbers[length(numbers)]

      sub_output$ExternalTarget = rep(target_id,nrow(sub_output))
      sub_output$ExternalTargetParameters = rep(external_target_parameters,nrow(sub_output))#rep(poorman::filter(theta_output,Target==target_id)$TargetParameters[1],nrow(sub_output))

      if (i==1)
      {
        lines_output = sub_output
      }
      else
      {
        lines_output = rbind(lines_output,sub_output)
      }
    }

    return(lines_output)
  }

  if (nesting_level==2)
  {
    stop("Currently function only written for nesting_level = 0 or 1.")
  }

  number_of_external_points = length(external_log_weights)
  if (as.mcmc && as.enk)
  {
    stop("Cannot treat output as both MCMC chains and an ensemble.")
  }

  if ( (ggmcmc==TRUE) && (as.mcmc==FALSE))
  {
    stop('as.mcmc must be set to TRUE if correct input for the ggmcmc package is to be produced.')
  }

  if ( (ggsmc==TRUE) && (as.mcmc==TRUE))
  {
    stop('as.mcmc must be set to FALSE if correct input for the ggsmc package is to be produced.')
  }

  if ( (ilike.output == TRUE) && (ggmcmc == TRUE) )
  {
    stop("ilike.output and ggmcmc cannot both be TRUE, since the output formats are incompatible.")
  }

  # ---- Read metadata from HDF5 output file ----
  h5_path <- file.path(results_directory, "output.h5")
  if (!file.exists(h5_path))
  {
    stop(paste("HDF5 output file not found:", h5_path))
  }
  h5f <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)

  # Root attributes: variable names and sizes
  variable_names <- h5f$attr_open("variable_names")$read()
  variable_sizes  <- as.numeric(h5f$attr_open("variable_sizes")$read())

  # Top-level datasets
  llhds_vec <- h5f[["log_likelihood"]]$read()
  times_vec  <- h5f[["time"]]$read()

  # output_lengths: HDF5 shape [n_iters x n_particles] (C row-major, h5dump convention).
  # hdf5r reverses dimension ordering: dims[1] = n_particles (C last dim),
  #                                    dims[2] = n_iters    (C first dim).
  ol_ds <- h5f[["output_lengths"]]
  ol_dims <- ol_ds$dims      # hdf5r: c(n_particles, n_iters)
  n_particles_dim <- ol_dims[1]
  n_iters <- ol_dims[2]
  ol_raw <- as.numeric(ol_ds$read())
  # hdf5r returns data in column-major / reversed order matching its dims.
  # Reshape to [n_iters x n_particles_dim].
  ol_mat <- matrix(ol_raw, nrow = n_iters, ncol = n_particles_dim, byrow = TRUE)
  output_lengths <- lapply(seq_len(n_iters), function(i) as.numeric(ol_mat[i, ]))

  # Convert to lists matching the old interface
  llhds <- as.list(llhds_vec)
  times  <- as.list(times_vec)

  if (length(variable_names) != length(variable_sizes))
  {
    stop("Variable names and sizes are of different lengths.")
  }

  output_names = rep("",sum(variable_sizes))
  output_index = rep("",sum(variable_sizes))
  counter = 1
  for (i in 1:length(variable_sizes))
  {
    for (j in 1:variable_sizes[i])
    {
      output_names[counter] = variable_names[i]
      output_index[counter] = j
      counter = counter + 1
    }
  }

  # Store the output in a data frame.
  all_dirs = order_directories(all_dirs)
  target_indices = get_directory_numbers(all_dirs)

  if (is.null(which.targets))
  {
    if (!as.mcmc)
    {
      which.targets = target_indices
    }
    else
    {
      which.targets = max(target_indices)
    }
  }
  else
  {
    inputted_which_targets = which.targets

    # target_indices = matrix(0,length(all_dirs))
    # for (k in 1:length(all_dirs))
    # {
    #   split = strsplit(all_dirs[k],"/")[[1]]
    #   target_indices[i] = strtoi(gsub("iteration","",split[length(split)]))
    # }

    which.targets = matrix(0,length(inputted_which_targets))
    for (k in 1:length(inputted_which_targets))
    {
      which.targets = which(target_indices==inputted_which_targets[k])
    }
  }

  for (k in which.targets)
  {
    #split = strsplit(all_dirs[k],"/")[[1]]
    #target = strtoi(gsub("iteration","",split[length(split)]))
    target = target_indices[k]

    iteration_directory = all_dirs[k]

    # Read vector_points from HDF5: stored as [n_rows x n_params] (C row-major).
    # hdf5r returns [n_params x n_rows] due to transposition, so we transpose back.
    iter_grp_path <- paste0("iteration/", target)
    vp_raw <- h5f[[iter_grp_path]][["vector_points"]]$read()
    if (is.matrix(vp_raw)) {
      output <- as.data.frame(t(vp_raw))
    } else {
      # hdf5r drops size-1 dimensions; vp_raw is a plain vector.
      # Reshape to [n_rows x n_cols] using the known number of parameters.
      output <- as.data.frame(matrix(vp_raw, ncol = sum(variable_sizes)))
    }

    old_column_names = names(output)
    number_of_points = nrow(output)
    number_of_importance_points = nrow(output) / number_of_external_points # number_of_importance_points in the waste-free SMC viewpoint
    # output_lengths[[target]] is the row for this iteration (1-indexed)
    chain_length = max(output_lengths[[target]])
    number_of_chains = number_of_importance_points/chain_length
    number_of_targets = length(all_dirs)

    output_names_column = rep(output_names,number_of_points)
    output_index_column = as.numeric(rep(output_index,number_of_points))

    # llhds and times are simple vectors in new HDF5 format (one entry per iteration)
    llhd_val = llhds[[target]]
    if (is.null(llhd_val)) llhd_val = 0
    llhd_column_matrix = sapply(lapply(1:number_of_external_points,FUN=function(i) { rep(llhd_val, ncol(output)*number_of_importance_points) } ),c)
    llhd_column = matrix(llhd_column_matrix,length(llhd_column_matrix))

    time_val = times[[target]]
    if (is.null(time_val)) time_val = 0
    time_column_matrix = sapply(lapply(1:number_of_external_points,FUN=function(i) { rep(time_val, ncol(output)*number_of_importance_points) } ),c)
    time_column = matrix(time_column_matrix,length(time_column_matrix))

    # number_of_points (=nrow(output)) equals number_of_external points multiplied by number of importance points

    schedule_parameters <- tryCatch(
      h5f[[iter_grp_path]][["schedule_parameters"]]$read(),
      error = function(e) NULL
    )
    if (is.null(schedule_parameters) || length(schedule_parameters) == 0)
    {
      schedule_parameters_column = matrix(0, nrow(output), 1)
      schedule_parameters <- NULL
    }
    else
    {
      schedule_parameters_column = matrix(as.character(schedule_parameters[1]), nrow(output) * ncol(output), 1)
    }

    # 1, ncol(output)*number_of_importance_points times
    # 2, ncol(output)*number_of_importance_points times
    # ...
    # number_of_external_points, ncol(output)*number_of_importance_points times
    external_column_matrix = sapply(lapply(1:number_of_external_points,FUN=function(i) { rep(i,ncol(output)*number_of_importance_points) } ),c)
    external_column = matrix(external_column_matrix,length(external_column_matrix))

    # all from this target
    target_column = rep(target,length(external_column))

    # (repeat each value (1:nchains) ncol(output)*niterations times)*number_of_external_points
    chain_fn = function(i) {rep(i,ncol(output)*chain_length)}
    chains_column = rep(sapply(lapply(1:number_of_chains,FUN=chain_fn),c),number_of_external_points)

    # (repeat each value (1:niterations) ncol(output) times)*number_of_external_points*nchains
    iteration_fn = function(i) {rep(i,ncol(output))}
    iterations_column = rep(sapply(lapply(1:chain_length,FUN=iteration_fn),c),number_of_external_points*number_of_chains)

    ess_val <- tryCatch(
      as.numeric(h5f[[iter_grp_path]][["ess"]]$read()),
      error = function(e) NA_real_
    )
    # ess is a single scalar per iteration; replicate to match long-format length.
    ess_column = rep(ess_val, length(iterations_column))

    output = poorman::pivot_longer(output,cols = poorman::everything(), values_to="Value")

    #external_column = matrix(0,nrow(output))
    #iterations_column = matrix(0,nrow(output),1)
    #chains_column = matrix(0,nrow(output),1)
    #target_column = matrix(0,nrow(output),1)

    # index = 1
    # for (p in 1:number_of_external_points)
    # {
    #   for (i in 1:length(output_lengths[[k]]))
    #   {
    #     for (j in 1:output_lengths[[k]][i])
    #     {
    #       #iterations_column[index] = j
    #       #chains_column[index] = i
    #       #target_column[index] = target
    #       #external_column[index] = p
    #       index = index + 1
    #     }
    #   }
    # }

    if (as.mcmc)
    {
      output = cbind(external_column,time_column,iterations_column,chains_column,output_names_column,output_index_column,output$Value)
      colnames(output) = c('ExternalIndex','Time','Iteration','Chain','ParameterName','Dimension',"Value")
    }
    else
    {
      if (as.enk)
      {
        if (is.null(schedule_parameters))
        {
          output = cbind(external_column,target_column,time_column,llhd_column,ess_column,iterations_column,chains_column,output_names_column,output_index_column,output$Value)
          colnames(output) = c('ExternalIndex','Target','Time','NormalisingConstant','ISESS','Iteration','Particle','ParameterName','Dimension',"Value")
        }
        else
        {
          output = cbind(external_column,target_column,time_column,llhd_column,ess_column,schedule_parameters_column,iterations_column,chains_column,output_names_column,output_index_column,output$Value)
          colnames(output) = c('ExternalIndex','Target','Time','NormalisingConstant','ISESS','TargetParameters','Iteration','Particle','ParameterName','Dimension',"Value")
        }
      }
      else
      {
        log_weight_raw <- h5f[[iter_grp_path]][["unnormalised_log_weights"]]$read()
        log_weight <- as.data.frame(matrix(as.numeric(log_weight_raw), ncol = 1))

        log_weights_for_each_external_point = lapply(1:number_of_external_points,FUN=function(i) { lw = nrow(log_weight)/number_of_external_points; o = matrix(0,lw); for (j in 1:lw) { o[j] = log_weight[j+lw*(i-1),] }; return(o); })
        normalised_weights_for_each_external_point = lapply(log_weights_for_each_external_point,FUN=function(i) { i-log_sum_exp(i) })
        normalised_weights_for_each_external_point = lapply(normalised_weights_for_each_external_point,FUN=function(i) { which_nan = which(is.nan(i)); o = i; o[which_nan] = -Inf; o })
        log_weight_list = lapply(1:number_of_external_points,FUN=function(i){ normalised_weights_for_each_external_point[[i]] + external_log_weights[i] })

        log_weight_fn = function(j) {matrix(sapply(1:length(log_weight_list[[j]]),FUN=function(i) { rep(log_weight_list[[j]][i],sum(variable_sizes)) }),sum(variable_sizes)*length(log_weight_list[[j]]))}

        log_weight_column = matrix(sapply(1:number_of_external_points,FUN=function(i) { log_weight_fn(i) }),length(iterations_column))

        #log_weight = log_weight$V1 - log_sum_exp(log_weight$V1)

        # (repeat each value of log_weight (1:nchains) ncol(output)*niterations times)*number_of_external_points
        #log_weight_fn = function(i) {rep(log_weight[i],sum(variable_sizes)*chain_length)}
        #log_weight_column = rep(sapply(lapply(1:number_of_chains,FUN=log_weight_fn),c),number_of_external_points)

        #external_log_weight_fn = function(i) {rep(external_log_weights[i],length(log_weight_column)/length(external_log_weights))}
        #external_log_weight_column = matrix(apply(lapply(1:number_of_external_points,FUN=log_weight_fn),c),length(log_weight_column))

        ancestor_index_raw <- tryCatch(
          as.integer(h5f[[iter_grp_path]][["ancestor_index"]]$read()) + 1L,
          error = function(e) seq_len(number_of_chains)
        )
        if (length(ancestor_index_raw) == 0) ancestor_index_raw <- seq_len(number_of_chains)
        ancestor_index <- matrix(ancestor_index_raw, ncol = 1)

        # (repeat each value of ancestor (1:nchains) ncol(output)*niterations times)*number_of_external_points
        ancestor_fn = function(i) {rep(ancestor_index[i,],sum(variable_sizes)*chain_length)}
        ancestor_index_column = rep(sapply(lapply(1:number_of_chains,FUN=ancestor_fn),c),number_of_external_points)

        if (is.null(schedule_parameters))
        {
          output = cbind(external_column,target_column,time_column,llhd_column,ess_column,iterations_column,chains_column,ancestor_index_column,log_weight_column,output_names_column,output_index_column,output$Value)
          colnames(output) = c('ExternalIndex','Target','Time','NormalisingConstant','ISESS','Iteration','Particle','AncestorIndex','LogWeight','ParameterName','Dimension',"Value")
        }
        else
        {
          output = cbind(external_column,target_column,time_column,llhd_column,ess_column,schedule_parameters_column,iterations_column,chains_column,ancestor_index_column,log_weight_column,output_names_column,output_index_column,output$Value)
          colnames(output) = c('ExternalIndex','Target','Time','NormalisingConstant','ISESS','TargetParameters','Iteration','Particle','AncestorIndex','LogWeight','ParameterName','Dimension',"Value")
        }

      }
    }

    output = as.data.frame(output)

    output$Dimension = as.integer(output$Dimension)
    output$ExternalIndex = as.integer(output$ExternalIndex)
    output$Value = as.numeric(output$Value)
    output$Iteration = as.integer(output$Iteration)

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

  if ("Target" %in% names(all_output))
  {
    all_output$Target = as.integer(all_output$Target)
  }

  if ("AncestorIndex" %in% names(all_output))
  {
    all_output$AncestorIndex = as.integer(all_output$AncestorIndex)
  }

  if ("LogWeight" %in% names(all_output))
  {
    all_output$LogWeight = as.numeric(all_output$LogWeight)
  }

  if ("Time" %in% names(all_output))
  {
    all_output$Time = as.numeric(all_output$Time)
  }

  if ("NormalisingConstant" %in% names(all_output))
  {
    all_output$NormalisingConstant = as.numeric(all_output$NormalisingConstant)
  }

  if ("ISESS" %in% names(all_output))
  {
    all_output$ISESS = as.numeric(all_output$ISESS)
  }

  if ("Particle" %in% names(all_output))
  {
    all_output$Particle = as.integer(all_output$Particle)
  }

  if ("Chain" %in% names(all_output))
  {
    all_output$Chain = as.integer(all_output$Chain)
  }

  if ( (ggmcmc==TRUE) && (as.mcmc==TRUE) )
  {
    new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },all_output$ParameterName,all_output$Dimension)
    all_output = subset(all_output,select = -c(ParameterName,Dimension))
    all_output$Parameter = new_variable_names

    attr(all_output,"nChains") = length(output_lengths[[k]])
    attr(all_output,"nParameters") = length(unique(new_variable_names))
    attr(all_output,"nIterations") = max(output_lengths[[k]])
    attr(all_output,"nBurnin") = 0
    attr(all_output,"nThin") = 1
    attr(all_output,"description") = description

    return(all_output)
  }

  if ( (ggsmc==FALSE) && (ilike.output==FALSE) )
  {
    all_output = ilike_pivot_wider(all_output)
  }

  return(all_output)
}

#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggmcmc (optional) Output in tidy format for plotting in ggmcmc package (default is FALSE if ilike.output is set to TRUE or not set and TRUE otherwise).
#' @param ilike.output (optional) Output can be processed by ilike,output package (default is TRUE if ggmcmc is set to FALSE or is not set and FALSE otherwise). Choosing TRUE for this argument is incompatible with ggmcmc=TRUE, since the two packages require different formatting of the output.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            ggmcmc = NULL,
                            ilike.output = NULL)
{
  if (is.null(ggmcmc) && is.null(ilike.output))
  {
    ggmcmc = FALSE
    ilike.output = TRUE
  }
  else if (!is.null(ggmcmc) && is.null(ilike.output))
  {
    if (ggmcmc == TRUE)
      ilike.output = FALSE
    else
      ilike.output = TRUE
  }
  else if (is.null(ggmcmc) && !is.null(ilike.output))
  {
    if (ilike.output == TRUE)
      ggmcmc = FALSE
    else
      ggmcmc = TRUE
  }
  else
  {
    if ( (ilike.output == TRUE) && (ggmcmc == TRUE) )
    {
      stop("ilike.output and ggmcmc cannot both be TRUE, since the output formats are incompatible.")
    }
  }

  output = load_smc_output(results_directory = results_directory,ggmcmc = ggmcmc,ggsmc = FALSE, ilike.output = ilike.output,as.mcmc=TRUE)
  if (ggmcmc==TRUE)
    names(output)[names(output) == 'Value'] = "value"
  return(output)
}

#' Loading ensemble Kalman output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggsmc (optional) Output in tidy format for plotting in ggsmc package.
#' @param external_log_weights (optional; for nested output only) The weights of the importance points external to the current folder. (default is 1, to be used at the top level of nested output)
#' @param external_target_parameters (optional; for nested output only) The parameters of the target external to the current folder. (default is "", corresponding to no parameters)
#' @param nesting_level (optional; for nested output only) The level of nesting at which to extract points. (default is 0, representing the top level of nested output)
#' @param factor (optional; for nested output only) The factor from which to extract points. (default is 0)
#' @return A list containing the ensemble members (called particles).
#' @export
load_enk_output = function(results_directory,
                           ggsmc = TRUE,
                           external_log_weights = c(0),
                           external_target_parameters = "",
                           nesting_level = 0,
                           factor = 1)
{
  return(load_smc_output(results_directory = results_directory,
                         ggsmc = ggsmc,
                         as.enk=TRUE,
                         external_log_weights = external_log_weights,
                         external_target_parameters = external_target_parameters,
                         nesting_level = nesting_level,
                         factor = factor))
}

#' #' Change data to be in "wide" format, with one column per variable.
#' #'
#' #' @param output The output from a "load_output" function.
#' #' @return The output in wide format.
#' ilike_pivot_wider = function(output)
#' {
#'   new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },output$ParameterName,output$Dimension)
#'   output = subset(output,select = -c(ParameterName,Dimension))
#'   output$Parameter = new_variable_names
#'   output = dplyr::distinct(output)
#'   return(tidyr::pivot_wider(output,names_from=Parameter,values_from=Value))
#' }
