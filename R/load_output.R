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
#' @param factor (optional; for nested output only) The factor from which to extract points. (default is 0)
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
                           factor = 0)
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
                                   external_log_weights = external_log_weights,#dplyr::filter(dplyr::filter(theta_output,Target==i),Dimension==1)$LogWeight,
                                   nesting_level = nesting_level-1)

      numbers = stringr::str_extract_all(all_dirs[i], "\\d+")
      target_id = as.integer(numbers[[1]][length(numbers[[1]])])
      sub_output$ExternalTarget = rep(target_id,nrow(sub_output))
      sub_output$ExternalTargetParameters = rep(external_target_parameters,nrow(sub_output))#rep(dplyr::filter(theta_output,Target==target_id)$TargetParameters[1],nrow(sub_output))

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

  # Make the output names.

  # New line for each iteration of SMC; one element in the line for each importance point -
  # tells us the length of the MCMC or single point. (since the file vector_points.txt is concatenated,
  # so we need this to split the file into different importance points). Same format for vector_variable_sizes and vector_variables.
  # When nested, a new line is written to these files for every external particle.

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

  llhds_file = file(paste(results_directory,"/log_likelihood.txt",sep=""),open="r")
  llhds = vector("list", length(variable_sizes))
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
      llhds[[counter]] = as.numeric(strsplit(line,split=" +")[[1]])
      #llhds[[counter]] = llhds[[counter]][!is.na(output_lengths[[counter]])]

      counter = counter + 1
    }
  }
  close(llhds_file)

  time_file = file(paste(results_directory,"/time.txt",sep=""),open="r")
  times = vector("list", length(variable_sizes))
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
      times[[counter]] = as.numeric(strsplit(line,split=" +")[[1]])
      #llhds[[counter]] = llhds[[counter]][!is.na(output_lengths[[counter]])]

      counter = counter + 1
    }
  }
  close(time_file)


  num_output_lengths_to_keep = length(output_lengths)/number_of_external_points
  output_lengths = output_lengths[1:num_output_lengths_to_keep]

  #number_of_points = length(output_lengths)

  if (length(variable_names)!=length(variable_sizes))
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

    output = utils::read.table(file=points_filename,header=FALSE)

    if (nrow(output)!=number_of_external_points*sum(output_lengths[[k]]))
    {
      stop("Number of rows in vector_points.txt file does not correspond to output_lengths.txt file.")
    }
    if (ncol(output)!=sum(variable_sizes))
    {
      stop("Number of columns in vector_points.txt file does not correspond to output_lengths.txt file.")
    }

    old_column_names = names(output)
    number_of_points = nrow(output)
    number_of_importance_points = nrow(output) / number_of_external_points # number_of_importance_points in the waste-free SMC viewpoint
    chain_length = max(output_lengths[[k]])
    number_of_chains = number_of_importance_points/chain_length

    output_names_column = rep(output_names,number_of_points)
    output_index_column = rep(output_index,number_of_points)

    llhd_column = rep(llhds[[k]],ncol(output)*number_of_importance_points*number_of_external_points)
    time_column = rep(times[[k]],ncol(output)*number_of_importance_points*number_of_external_points)

    # number_of_points (=nrow(output)) equals number_of_external points multiplied by number of importance points

    schedule_parameters_filename = paste(iteration_directory,"/schedule_parameters.txt",sep="")
    tryCatch( {schedule_parameters = utils::read.table(file=schedule_parameters_filename,header=FALSE) }
              , error = function(e) {schedule_parameters <<- NULL})
    if (is.null(schedule_parameters))
    {
      schedule_parameters_column = matrix(0,nrow(output),1)
    }
    else
    {
      schedule_parameters_column = matrix(paste(schedule_parameters[1,]),nrow(output)*ncol(output),1)
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

    ess_filename = paste(iteration_directory,"/ess.txt",sep="")
    ess = utils::read.table(file=ess_filename,header=FALSE,sep=",")

    ess_for_each_external_point = lapply(1:number_of_external_points,FUN=function(i) { lw = nrow(ess)/number_of_external_points; o = matrix(0,lw); for (j in 1:lw) { o[j] = ess[j+lw*(i-1),] }; return(o); })
    ess_list = lapply(1:number_of_external_points,FUN=function(i){ ess_for_each_external_point[[i]] })
    ess_fn = function(j) {matrix(sapply(1:length(ess_list[[j]]),FUN=function(i) { rep(ess_list[[j]][i],sum(variable_sizes)) }),sum(variable_sizes)*length(ess_list[[j]]))}
    ess_column = matrix(sapply(1:number_of_external_points,FUN=function(i) { ess_fn(i) }),length(iterations_column))

    output = tidyr::pivot_longer(output,cols = dplyr::everything(), values_to="Value")

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
        log_weight_filename = paste(iteration_directory,"/unnormalised_log_weights.txt",sep="")
        log_weight = utils::read.table(file=log_weight_filename,header=FALSE,sep=",")

        log_weights_for_each_external_point = lapply(1:number_of_external_points,FUN=function(i) { lw = nrow(log_weight)/number_of_external_points; o = matrix(0,lw); for (j in 1:lw) { o[j] = log_weight[j+lw*(i-1),] }; return(o); })
        normalised_weights_for_each_external_point = lapply(log_weights_for_each_external_point,FUN=function(i) { i-log_sum_exp(i) })
        normalised_weights_for_each_external_point = lapply(normalised_weights_for_each_external_point,FUN=function(i) { which_nan = which(is.nan(i)); o = i; o[which_nan] = -Inf; o })
        log_weight_list = lapply(1:number_of_external_points,FUN=function(i){ normalised_weights_for_each_external_point[[i]] + external_log_weights[i] })

        log_weight_fn = function(j) {matrix(sapply(1:length(log_weight_list[[j]]),FUN=function(i) { rep(log_weight_list[[j]][i],sum(variable_sizes)) }),sum(variable_sizes)*length(log_weight_list[[j]]))}

        log_weight_column = matrix(sapply(1:number_of_external_points,FUN=function(i) { log_weight_fn(i) }),length(iterations_column))

        #browser()
        #log_weight = log_weight$V1 - log_sum_exp(log_weight$V1)

        # (repeat each value of log_weight (1:nchains) ncol(output)*niterations times)*number_of_external_points
        #log_weight_fn = function(i) {rep(log_weight[i],sum(variable_sizes)*chain_length)}
        #log_weight_column = rep(sapply(lapply(1:number_of_chains,FUN=log_weight_fn),c),number_of_external_points)

        #external_log_weight_fn = function(i) {rep(external_log_weights[i],length(log_weight_column)/length(external_log_weights))}
        #external_log_weight_column = matrix(apply(lapply(1:number_of_external_points,FUN=log_weight_fn),c),length(log_weight_column))

        ancestor_index_filename = paste(iteration_directory,"/ancestor_index.txt",sep="")
        tryCatch( {ancestor_index = utils::read.table(file=ancestor_index_filename,header=FALSE,sep=",") + 1 }
                  , error = function(e) {ancestor_index <<- matrix(1:number_of_chains,number_of_chains)})

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
    new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },output$ParameterName,output$Dimension)
    output = subset(output,select = -c(ParameterName,Dimension))
    output$Parameter = new_variable_names

    attr(output,"nChains") = length(output_lengths[[k]])
    attr(output,"nParameters") = length(unique(new_variable_names))
    attr(output,"nIterations") = max(output_lengths[[k]])
    attr(output,"nBurnin") = 0
    attr(output,"nThin") = 1
    attr(output,"description") = description

    return(output)
  }

  if ( (ggsmc==FALSE) && (ilike.output==FALSE) )
  {
    new_variable_names = mapply(FUN = function(a,b) { paste(a,"_",b,sep="") },all_output$ParameterName,all_output$Dimension)
    all_output = subset(all_output,select = -c(ParameterName,Dimension))
    all_output$Parameter = new_variable_names
    output = dplyr::distinct(output)
    all_output = tidyr::pivot_wider(all_output, names_from = "Parameter", values_from = "Value")
  }

  return(all_output)
}

#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggmcmc (optional) Output in tidy format for plotting in ggmcmc package (default is TRUE).
#' @param ilike.output (optional) Output can be processed by ilike,output package (default is FALSE). Choosing TRUE for this argument is incompatiable with ggmcmc=TRUE, since the two packages require different formatting of the output.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            ggmcmc = NULL,
                            ilike.output = NULL)
{
  if (is.null(ggmcmc) && is.null(ilike.output))
  {
    ggmcmc = TRUE
    ilike.output = FALSE
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

  return(load_smc_output(results_directory = results_directory,ggmcmc = ggmcmc,ggsmc = FALSE, ilike.output = ilike.output,as.mcmc=TRUE))
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
                           factor = 0)
{
  return(load_smc_output(results_directory = results_directory,
                         ggsmc = ggsmc,
                         as.enk=TRUE,
                         external_log_weights = external_log_weights,
                         external_target_parameters = external_target_parameters,
                         nesting_level = nesting_level,
                         factor = factor))
}
