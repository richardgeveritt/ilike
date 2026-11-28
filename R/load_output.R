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

#' Change ilike output from the wide format (one row per point, one column per
#' scalar, as given by \code{wide = TRUE} in the load functions or by
#' \code{ilike_pivot_wider}) to the long format (one row per scalar, with
#' columns ParameterName, Dimension and Value), as used by ggsmc.
#'
#' @param output Output in the wide format.
#' @return The output in the long format: the index columns, then ParameterName, Dimension and Value; rows ordered point by point, as in the load functions.
#' @export
ilike_pivot_longer <- function(output)
{
  index_names = c("ExternalIndex", "ExternalTarget", "ExternalTargetParameters", "Target", "Time",
                  "NormalisingConstant", "ISESS", "TargetParameters", "Iteration", "Particle",
                  "Chain", "AncestorIndex", "LogWeight", "Sign")
  output = as.data.frame(output)
  is_index = names(output) %in% index_names
  value_names = names(output)[!is_index]
  if (length(value_names) == 0)
    stop("ilike_pivot_longer: no value columns (named variable_dimension, e.g. theta_1) found.")
  parts = regmatches(value_names, regexec("^(.*)_([0-9]+)$", value_names))
  if (any(lengths(parts) != 3))
    stop("ilike_pivot_longer: columns not named variable_dimension: ",
         paste(value_names[lengths(parts) != 3], collapse = ", "))
  parameter_names = vapply(parts, function(x) x[2], character(1))
  dimensions = as.integer(vapply(parts, function(x) x[3], character(1)))
  n = nrow(output); p = length(value_names)
  long = output[rep(seq_len(n), each = p), is_index, drop = FALSE]
  rownames(long) = NULL
  long$ParameterName = rep(parameter_names, times = n)
  long$Dimension = rep(dimensions, times = n)
  long$Value = as.vector(t(as.matrix(output[, value_names, drop = FALSE])))
  long
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
#' @param wide (optional) If TRUE, return one row per point, with one column per scalar (named variable_dimension, e.g. theta_1) next to the index columns, instead of one row per scalar (default FALSE). Much faster for large outputs; ilike_pivot_longer converts it to the long format (needed by ggsmc).
#' @return A data frame containing the SMC output.
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
                           factor = 1,
                           wide = FALSE)
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
  # The iterationK directories hold the factor outputs; they are absent when
  # the algorithm was run with write_factor_output = FALSE.
  if (length(all_dirs) > 0)
    all_dirs = order_directories(all_dirs)
  else if (nesting_level > 0)
    stop("No factor outputs in ", results_directory, ": was the algorithm run with write_factor_output = FALSE?")

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
    # The outer points and their weights: each outer iteration's points, in
    # the order their nested outputs were written (point by point), give the
    # external weights of that iteration's nested output.
    outer = load_smc_output(results_directory = description, as.mcmc = as.mcmc,
                            ggsmc = !as.mcmc, directory_prefix = directory_prefix, wide = TRUE)
    lines_output = NULL
    for (i in 1:length(all_dirs))
    {
      numbers = extract_numbers(all_dirs[i])
      target_id = numbers[length(numbers)]
      if (!dir.exists(file.path(all_dirs[i], paste(factor_name, "_smc", sep = ""))))
        next
      # MCMC output: the outer points loaded are those of the final iteration
      if (as.mcmc && (target_id != max(get_directory_numbers(all_dirs))))
        next
      if ("Target" %in% names(outer))
        outer_points = outer[outer$Target == target_id, , drop = FALSE]
      else
        outer_points = outer
      if ("LogWeight" %in% names(outer_points))
        external_log_weights = outer_points$LogWeight
      else
        external_log_weights = rep(0, nrow(outer_points))
      if ("TargetParameters" %in% names(outer_points))
        external_target_parameters = as.character(outer_points$TargetParameters[1])
      sub_output = load_smc_output(results_directory=all_dirs[i],
                                   ggmcmc=FALSE,
                                   ggsmc=ggsmc,
                                   as.mcmc = FALSE,
                                   as.enk = as.enk,
                                   which.targets = which.targets,
                                   directory_prefix=factor_name,
                                   external_log_weights = external_log_weights,#poorman::filter(poorman::filter(theta_output,Target==i),Dimension==1)$LogWeight,
                                   nesting_level = nesting_level-1,
                                   wide = wide)


      sub_output$ExternalTarget = rep(target_id,nrow(sub_output))
      sub_output$ExternalTargetParameters = rep(external_target_parameters,nrow(sub_output))#rep(poorman::filter(theta_output,Target==target_id)$TargetParameters[1],nrow(sub_output))

      if (is.null(lines_output))
        lines_output = sub_output
      else
        lines_output = rbind(lines_output,sub_output)
    }

    if (is.null(lines_output))
      stop("No nested output for factor ", factor, " in ", results_directory, ".")
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

  if (wide && ggmcmc)
    stop("wide = TRUE and ggmcmc = TRUE are incompatible: ggmcmc needs the long format.")
  # ggsmc = FALSE and ilike.output = FALSE asked for the wide format
  # (previously made by ilike_pivot_wider at the end).
  if ( (ggsmc==FALSE) && (ilike.output==FALSE) && (ggmcmc==FALSE) )
    wide = TRUE

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
  # A nested inner file (written by SMCOutput::write_compact_nested) holds
  # the outputs of several outer points concatenated: output lengths per
  # iteration group, and one log-likelihood, time and ESS per outer point.
  nested_file <- "nested_external_points" %in% hdf5r::h5attr_names(h5f)
  output_lengths <- list()
  if (!nested_file)
  {
    ol_ds <- h5f[["output_lengths"]]
    ol_dims <- ol_ds$dims      # hdf5r: c(n_particles, n_iters)
    n_particles_dim <- ol_dims[1]
    n_iters <- ol_dims[2]
    ol_raw <- as.numeric(ol_ds$read())
    # hdf5r returns data in column-major / reversed order matching its dims.
    # Reshape to [n_iters x n_particles_dim].
    ol_mat <- matrix(ol_raw, nrow = n_iters, ncol = n_particles_dim, byrow = TRUE)
    output_lengths <- lapply(seq_len(n_iters), function(i) as.numeric(ol_mat[i, ]))
  }

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

  # Store the output in a data frame. The stored iterations are the groups
  # iteration/K of output.h5 (the matching iterationK directories exist only
  # when factor outputs were written).
  target_indices = sort(as.integer(names(h5f[["iteration"]])))
  all_dirs = file.path(results_directory, paste0("iteration", target_indices))

  # which.targets holds positions in target_indices (the loop below reads
  # target_indices[k]). It used to hold the target numbers themselves, which
  # only worked while the stored iterations were exactly 1, 2, ...; and a
  # given which.targets kept only its last target.
  if (is.null(which.targets))
  {
    if (!as.mcmc)
      which.targets = seq_along(target_indices)
    else
      which.targets = which.max(target_indices)
  }
  else
  {
    which.targets = which(target_indices %in% which.targets)
    if (length(which.targets) == 0)
      stop("load_smc_output: none of the requested targets is in the output.")
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
    if (nested_file)
      chain_length = max(as.numeric(h5f[[iter_grp_path]][["output_lengths"]]$read()))
    else
      chain_length = max(output_lengths[[target]])
    number_of_chains = number_of_importance_points/chain_length
    number_of_targets = length(all_dirs)

    if (nested_file)
    {
      # one value per outer (external) point
      llhd_val = as.numeric(llhds_vec)
      time_val = as.numeric(times_vec)
    }
    else
    {
      llhd_val = llhds[[target]]
      if (is.null(llhd_val)) llhd_val = 0
      time_val = times[[target]]
      if (is.null(time_val)) time_val = 0
    }
    schedule_parameters <- tryCatch(
      h5f[[iter_grp_path]][["schedule_parameters"]]$read(),
      error = function(e) NULL
    )
    if (is.null(schedule_parameters) || length(schedule_parameters) == 0)
      schedule_parameters <- NULL

    # Build the index columns once per point (one row of `output` per
    # point), then either keep one row per point (wide) or repeat each index
    # value once per scalar (long). Rows are point-major: external point,
    # then chain (particle), then stored iteration. (This replaces a cbind
    # of mixed types, which went through a character matrix: slow, and
    # numbers kept only 15 significant digits.)
    p_total = sum(variable_sizes)
    points_per_external = number_of_chains*chain_length
    point_external = rep(seq_len(number_of_external_points), each = points_per_external)
    point_chain = rep(rep(seq_len(number_of_chains), each = chain_length), times = number_of_external_points)
    point_iteration = rep(seq_len(chain_length), times = number_of_chains*number_of_external_points)
    llhd_val = as.numeric(llhd_val)
    time_val = as.numeric(time_val)
    if (nested_file)
    {
      if ((length(llhd_val) != number_of_external_points) || (length(time_val) != number_of_external_points))
        stop("load_smc_output: the nested output holds ", length(llhd_val), " outer points, but ",
             number_of_external_points, " external weights were given.")
      llhd_val = llhd_val[point_external]
      time_val = time_val[point_external]
    }

    ess_val <- tryCatch(
      as.numeric(h5f[[iter_grp_path]][["ess"]]$read()),
      error = function(e) NA_real_
    )
    if (nested_file && (length(ess_val) == number_of_external_points))
      ess_val = ess_val[point_external]

    index = data.frame(ExternalIndex = point_external)
    if (as.mcmc)
    {
      index$Time = time_val
      index$Iteration = point_iteration
      index$Chain = point_chain
      # (dev/plan.md 3.22 U9 R2) signed pseudo-marginal MCMC: each point's sign
      if ("signs" %in% names(h5f[[iter_grp_path]]))
      {
        sign_raw = as.numeric(h5f[[iter_grp_path]][["signs"]]$read())
        if (length(sign_raw) != number_of_points)
          stop("load_mcmc_output: the output holds ", length(sign_raw), " signs for ", number_of_points, " points.")
        index$Sign = sign_raw
      }
    }
    else
    {
      index$Target = rep(as.integer(target), number_of_points)
      # (dev/plan.md 3.22 U9 R5) IS and SMC with factor estimates that may be
      # negative: each point's sign
      if ("signs" %in% names(h5f[[iter_grp_path]]))
      {
        sign_raw = as.numeric(h5f[[iter_grp_path]][["signs"]]$read())
        if (length(sign_raw) != number_of_points)
          stop("load_smc_output: the output holds ", length(sign_raw), " signs for ", number_of_points, " points.")
        index$Sign = sign_raw
      }
      index$Time = time_val
      index$NormalisingConstant = llhd_val
      index$ISESS = ess_val
      if (!is.null(schedule_parameters))
        index$TargetParameters = as.character(schedule_parameters[1])
      index$Iteration = point_iteration
      index$Particle = point_chain
      if (!as.enk)
      {
        log_weight_raw <- as.numeric(h5f[[iter_grp_path]][["unnormalised_log_weights"]]$read())
        lw_per_external = length(log_weight_raw)/number_of_external_points
        log_weight_list = lapply(seq_len(number_of_external_points), FUN = function(i) {
          lw = log_weight_raw[(i - 1)*lw_per_external + seq_len(lw_per_external)]
          lw = lw - log_sum_exp(lw)
          lw[is.nan(lw)] = -Inf
          lw + external_log_weights[i]
        })
        # each particle's weight covers its chain_length stored iterations
        point_log_weight = unlist(lapply(log_weight_list, FUN = function(w) rep(as.numeric(w), each = chain_length)))
        if (length(point_log_weight) != number_of_points)
        {
          stop(paste("load_smc_output: target ", target, " has ", length(point_log_weight)/chain_length,
                     " weights but ", number_of_chains*number_of_external_points, " particles.", sep = ""))
        }

        ancestor_index_raw <- tryCatch(
          as.integer(h5f[[iter_grp_path]][["ancestor_index"]]$read()) + 1L,
          error = function(e) seq_len(number_of_chains)
        )
        if (length(ancestor_index_raw) == 0) ancestor_index_raw <- seq_len(number_of_chains)
        if (length(ancestor_index_raw) == number_of_chains*number_of_external_points)
          # one set of ancestors per external point (nested output)
          index$AncestorIndex = rep(ancestor_index_raw, each = chain_length)
        else
          index$AncestorIndex = rep(rep(ancestor_index_raw[seq_len(number_of_chains)], each = chain_length), times = number_of_external_points)
        index$LogWeight = point_log_weight
      }
    }

    values = as.matrix(output)
    dimnames(values) = NULL
    if (wide)
    {
      colnames(values) = paste(output_names, output_index, sep = "_")
      output = cbind(index, as.data.frame(values, stringsAsFactors = FALSE))
    }
    else
    {
      output = index[rep(seq_len(number_of_points), each = p_total), , drop = FALSE]
      rownames(output) = NULL
      output$ParameterName = rep(output_names, times = number_of_points)
      output$Dimension = as.integer(rep(output_index, times = number_of_points))
      output$Value = as.vector(t(values))
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

  return(all_output)
}

#' Loading MCMC output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggmcmc (optional) Output in tidy format for plotting in ggmcmc package (default is FALSE if ilike.output is set to TRUE or not set and TRUE otherwise).
#' @param ilike.output (optional) Output can be processed by ilike,output package (default is TRUE if ggmcmc is set to FALSE or is not set and FALSE otherwise). Choosing TRUE for this argument is incompatible with ggmcmc=TRUE, since the two packages require different formatting of the output.
#' @param wide (optional) If TRUE, return one row per point, with one column per scalar (named variable_dimension, e.g. theta_1), instead of one row per scalar (default FALSE). Much faster for large outputs; ilike_pivot_longer converts it to the long format.
#' @return A list containing the MCMC chains.
#' @export
load_mcmc_output = function(results_directory,
                            ggmcmc = NULL,
                            ilike.output = NULL,
                            wide = FALSE)
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

  output = load_smc_output(results_directory = results_directory,ggmcmc = ggmcmc,ggsmc = FALSE, ilike.output = ilike.output,as.mcmc=TRUE,wide=wide)
  output = mask_padded_entries(output, results_directory, wide = wide, ggmcmc = ggmcmc)
  if (ggmcmc==TRUE)
    names(output)[names(output) == 'Value'] = "value"
  return(output)
}

# (dev/plan.md 3.22 U6, C314) A model with an indicator k ~ categorical(…)
# and variables whose loops run to k stores those variables padded to k's
# largest value. MCMC records which they are in the output file, and
# load_mcmc_output marks the entries after k's value at each point as missing.
padded_output_file <- function(results_directory)
{
  # (the results folder, or, inside MCMC, the prefix the writer appends "_smc" to)
  for (path in c(file.path(results_directory, "ilike_smc", "output.h5"), file.path(paste0(results_directory, "_smc"), "output.h5"),
                 file.path(results_directory, "output.h5")))
    if (file.exists(path))
      return(path)
  ""
}

write_padded_variables <- function(recipe, results_directory)
{
  if ((results_directory == "") || (length(recipe$padded_variables) == 0))
    return(invisible(NULL))
  path = padded_output_file(results_directory)
  if (path == "")
    return(invisible(NULL))
  h5f = hdf5r::H5File$new(path, mode = "r+")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)
  hdf5r::h5attr(h5f, "padded_variables") = as.character(recipe$padded_variables)
  hdf5r::h5attr(h5f, "padded_indicators") = as.character(recipe$padded_indicators)
  invisible(NULL)
}

mask_padded_entries <- function(output, results_directory, wide = FALSE, ggmcmc = FALSE)
{
  path = padded_output_file(results_directory)
  if (path == "")
    return(output)
  h5f = hdf5r::H5File$new(path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)
  if (!("padded_variables" %in% hdf5r::h5attr_names(h5f)))
    return(output)
  variables = h5f$attr_open("padded_variables")$read()
  indicators = h5f$attr_open("padded_indicators")$read()
  keys = intersect(c("ExternalIndex", "Chain", "Iteration"), names(output))
  key_of = function(rows) do.call(paste, c(lapply(keys, function(k) output[[k]][rows]), sep = "_"))
  for (j in seq_along(variables))
  {
    v = variables[j]; k = indicators[j]
    if (wide)
    {
      k_column = intersect(c(paste0(k, "_1"), k), names(output))
      if (length(k_column) == 0)
        next
      kk = output[[k_column[1]]]
      columns = grep(paste0("^", v, "_[0-9]+$"), names(output), value = TRUE)
      for (column in columns)
      {
        d = as.numeric(sub(paste0("^", v, "_"), "", column))
        output[[column]][d > kk] = NA
      }
      next
    }
    if (ggmcmc)
    {
      k_rows = which(output$Parameter == paste0(k, "_1"))
      v_rows = grep(paste0("^", v, "_[0-9]+$"), output$Parameter)
      d = as.numeric(sub(paste0("^", v, "_"), "", output$Parameter[v_rows]))
    }
    else
    {
      k_rows = which(output$ParameterName == k)
      v_rows = which(output$ParameterName == v)
      d = output$Dimension[v_rows]
    }
    if ((length(k_rows) == 0) || (length(v_rows) == 0))
      next
    kk = output$Value[k_rows][match(key_of(v_rows), key_of(k_rows))]
    inactive = !is.na(kk) & (d > kk)
    output$Value[v_rows[inactive]] = NA
  }
  output
}

#' Loading ensemble Kalman output into R memory.
#'
#' @param results_directory The folder in which the results are stored.
#' @param ggsmc (optional) Output in tidy format for plotting in ggsmc package.
#' @param external_log_weights (optional; for nested output only) The weights of the importance points external to the current folder. (default is 1, to be used at the top level of nested output)
#' @param external_target_parameters (optional; for nested output only) The parameters of the target external to the current folder. (default is "", corresponding to no parameters)
#' @param nesting_level (optional; for nested output only) The level of nesting at which to extract points. (default is 0, representing the top level of nested output)
#' @param factor (optional; for nested output only) The factor from which to extract points. (default is 0)
#' @param wide (optional) If TRUE, return one row per point, with one column per scalar (named variable_dimension, e.g. theta_1), instead of one row per scalar (default FALSE). Much faster for large outputs; ilike_pivot_longer converts it to the long format.
#' @return A list containing the ensemble members (called particles).
#' @export
load_enk_output = function(results_directory,
                           ggsmc = TRUE,
                           external_log_weights = c(0),
                           external_target_parameters = "",
                           nesting_level = 0,
                           factor = 1,
                           wide = FALSE)
{
  return(load_smc_output(results_directory = results_directory,
                         ggsmc = ggsmc,
                         as.enk=TRUE,
                         wide = wide,
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

#' Loading the adaptation of adaptive MCMC moves into R memory.
#'
#' For MCMC moves with an adaptor (an m_adaptor or mh_adaptor block, or an
#' adapt_m_proposal or adapt_mh_proposal function), each chain records the
#' adapted parameters of its proposal: at every stored iteration, the scales
#' and the running estimate of the target's mean (the trace); and at the end
#' of the chain, also the estimated covariance (its upper triangle) and, for
#' principal-components moves, the eigenvalues and directions (the final
#' state). Parameter names are prefixed with "moveK." when the MCMC consists
#' of several moves.
#'
#' @param results_directory The folder in which the results are stored.
#' @param target (optional) The stored SMC iteration to read (default: the last one, which for MCMC holds the chains).
#' @return A list with data frames trace (columns Chain, Iteration, Parameter, Value; Iteration counts stored iterations, as in load_mcmc_output) and final (columns Chain, Parameter, Value).
#' @export
load_mcmc_adaptation = function(results_directory,
                                target = NULL)
{
  h5_path <- file.path(results_directory, "ilike_smc", "output.h5")
  if (!file.exists(h5_path))
    stop(paste("HDF5 output file not found:", h5_path))
  h5f <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)
  if (!("adaptation_trace_names" %in% hdf5r::h5attr_names(h5f)))
    stop("load_mcmc_adaptation: the output holds no adaptation (no adaptive MCMC move was used).")
  trace_names <- h5f$attr_open("adaptation_trace_names")$read()
  final_names <- h5f$attr_open("adaptation_final_names")$read()

  targets <- sort(as.integer(names(h5f[["iteration"]])))
  if (is.null(target))
    target <- max(targets)
  group <- h5f[[paste0("iteration/", target)]]
  if (!("adaptation_trace" %in% names(group)))
    stop("load_mcmc_adaptation: iteration ", target, " holds no adaptation.")

  # stored as [rows x columns]; hdf5r returns the transpose (and drops
  # dimensions of size 1)
  read_rows <- function(dataset, number_of_columns)
  {
    raw <- dataset$read()
    if (is.matrix(raw))
      return(t(raw))
    matrix(raw, ncol = number_of_columns, byrow = TRUE)
  }
  trace <- read_rows(group[["adaptation_trace"]], length(trace_names))
  final <- read_rows(group[["adaptation_final"]], length(final_names))

  # the stored points of chain k are rows sum(lengths[1:(k-1)]) + 1, ...
  ol <- h5f[["output_lengths"]]
  ol_dims <- ol$dims
  ol_mat <- matrix(as.numeric(ol$read()), nrow = ol_dims[2], ncol = ol_dims[1], byrow = TRUE)
  lengths <- ol_mat[which(targets == target), ]
  if (sum(lengths) != nrow(trace))
    stop("load_mcmc_adaptation: the adaptation trace does not match the stored points.")
  chain <- rep(seq_along(lengths), lengths)
  iteration <- sequence(lengths)

  list(trace = data.frame(Chain = rep(chain, times = length(trace_names)),
                          Iteration = rep(iteration, times = length(trace_names)),
                          Parameter = rep(trace_names, each = nrow(trace)),
                          Value = as.vector(trace),
                          stringsAsFactors = FALSE),
       final = data.frame(Chain = rep(seq_len(nrow(final)), times = length(final_names)),
                          Parameter = rep(final_names, each = nrow(final)),
                          Value = as.vector(final),
                          stringsAsFactors = FALSE))
}


#' Load the tuning values of an SMC run with esjd adaptation.
#'
#' With `smc_adapt = esjd(...)` on a move, each particle of `SMC_with_MCMC` moves
#' with its own tuning values, chosen from how far the particles moved at the
#' previous step. This function reads the summary of those values written at
#' each step: per tuned move and parameter, the mean, standard deviation and
#' 10%, 50% and 90% quantiles over the particles.
#'
#' @param results_directory The results directory of the run (`results_name`).
#' @return A data frame with one row per move step and tuner: columns step (the move step: 1 for the first move),
#'   tuner (its number, when several moves or proposals are tuned), variable, and for each tuned parameter
#'   (scale, or step and steps) its mean, sd, q10, q50 and q90; for pretune also step_max and steps_max (the
#'   ranges of the next trial moves); for choose, kernelK.share (the share of the particles using each proposal).
#' @export
load_smc_tuning = function(results_directory)
{
  h5_path = file.path(results_directory, "ilike_smc", "output.h5")
  h5f = if (file.exists(h5_path)) hdf5r::H5File$new(h5_path, mode = "r") else NULL
  if (is.null(h5f) || !h5f$exists("tuning"))
  {
    if (!is.null(h5f)) h5f$close_all()
    stop("load_smc_tuning: no tuning values in ", results_directory, " (the run used no esjd adaptation).")
  }
  on.exit(h5f$close_all(), add = TRUE)
  # one group per tuner: the variable as an attribute, the step and each
  # statistic as datasets with one value per move step
  tuners = h5f[["tuning"]]$ls()$name
  tuners = tuners[order(as.integer(tuners))]
  parts = lapply(tuners, function(k)
  {
    g = h5f[["tuning"]][[k]]
    statistics = g$attr_open("statistics")$read()
    step = as.numeric(g[["step"]]$read())
    frame = data.frame(step = step, tuner = as.integer(k), variable = g$attr_open("variable")$read(),
                       stringsAsFactors = FALSE)
    for (s in statistics)
      frame[[s]] = as.numeric(g[[s]]$read())
    frame
  })
  # one row per step and tuner, one column per statistic (NA where a tuner has none)
  statistics = unique(unlist(lapply(parts, function(f) setdiff(names(f), c("step", "tuner", "variable")))))
  parts = lapply(parts, function(f) { for (s in setdiff(statistics, names(f))) f[[s]] = NA_real_; f[, c("step", "tuner", "variable", statistics)] })
  output = do.call(rbind, parts)
  output = output[order(output$tuner, output$step), ]
  rownames(output) = NULL
  output
}

#' Load the number of particles of SMC^2's filters at each step.
#'
#' In SMC^2 (an outer `SMC_with_MCMC` whose `sequence = data(...)` takes in the
#' observations of a factor `pf(inner, ...)`), the setting
#' `increase = exchange(below = b, multiply = m, max = M)` of `pf` gives the
#' filters more particles when the acceptance rate of the moves falls below
#' `b` (the exchange step of Chopin, Jacob and Papaspiliopoulos, 2013). This
#' function reads the record of each step that moved the particles.
#'
#' @param results_directory The results directory of the run (`results_name`).
#' @return A data frame with one row per step that moved the particles and filter with `exchange`: columns step,
#'   factor (the factor's number in the model), acceptance_rate (the fraction of particles the step's moves
#'   changed), particles_before and particles_after (the filters' number of particles before and after the step's
#'   exchange; equal when there was none).
#' @export
load_smc_inner_particles = function(results_directory)
{
  h5_path = file.path(results_directory, "ilike_smc", "output.h5")
  h5f = if (file.exists(h5_path)) hdf5r::H5File$new(h5_path, mode = "r") else NULL
  if (is.null(h5f) || !h5f$exists("inner_particles"))
  {
    if (!is.null(h5f)) h5f$close_all()
    stop("load_smc_inner_particles: no record in ", results_directory, " (no filter with increase = exchange(...)).")
  }
  on.exit(h5f$close_all(), add = TRUE)
  factors = h5f[["inner_particles"]]$ls()$name
  factors = factors[order(as.integer(factors))]
  output = do.call(rbind, lapply(factors, function(k)
  {
    g = h5f[["inner_particles"]][[k]]
    data.frame(step = as.numeric(g[["step"]]$read()), factor = as.integer(k),
               acceptance_rate = as.numeric(g[["acceptance_rate"]]$read()),
               particles_before = as.numeric(g[["particles_before"]]$read()),
               particles_after = as.numeric(g[["particles_after"]]$read()))
  }))
  output = output[order(output$factor, output$step), ]
  rownames(output) = NULL
  output
}

#' Load the per-iteration diagnostics of an MCMC run with nuts.
#'
#' A `nuts` move records, at every stored iteration of each chain, the depth
#' of the tree it built (the trajectory doubled that many times), the number
#' of leapfrog steps (gradient evaluations), whether the trajectory diverged
#' (the error in the Hamiltonian passed 1000, a sign that the step is too
#' large for the curvature of the target there), and the acceptance statistic
#' (the average acceptance probability over the last doubling, which the
#' adaptors drive towards their target). Many divergences after warm-up mean
#' the chain may miss part of the target; a tree depth often at max_depth
#' means the step is too small or the target poorly scaled.
#'
#' @param results_directory The folder in which the results are stored.
#' @param target (optional) The stored SMC iteration to read (default: the last one, which for MCMC holds the chains).
#' @return A data frame with one row per stored iteration: Chain, Iteration (counting stored iterations, as in load_mcmc_output), tree_depth, leapfrog_steps, divergent (logical) and accept_stat.
#' @export
load_mcmc_diagnostics = function(results_directory,
                                 target = NULL)
{
  h5_path <- file.path(results_directory, "ilike_smc", "output.h5")
  if (!file.exists(h5_path))
    stop(paste("HDF5 output file not found:", h5_path))
  h5f <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(try(h5f$close_all(), silent = TRUE), add = TRUE)
  if (!("diagnostics_names" %in% hdf5r::h5attr_names(h5f)))
    stop("load_mcmc_diagnostics: the output holds no diagnostics (no nuts move was used).")
  names <- h5f$attr_open("diagnostics_names")$read()

  targets <- sort(as.integer(names(h5f[["iteration"]])))
  if (is.null(target))
    target <- max(targets)
  group <- h5f[[paste0("iteration/", target)]]
  if (!("diagnostics" %in% names(group)))
    stop("load_mcmc_diagnostics: iteration ", target, " holds no diagnostics.")
  raw <- group[["diagnostics"]]$read()
  values <- if (is.matrix(raw)) t(raw) else matrix(raw, ncol = length(names), byrow = TRUE)

  ol <- h5f[["output_lengths"]]
  ol_dims <- ol$dims
  ol_mat <- matrix(as.numeric(ol$read()), nrow = ol_dims[2], ncol = ol_dims[1], byrow = TRUE)
  lengths <- ol_mat[which(targets == target), ]
  if (sum(lengths) != nrow(values))
    stop("load_mcmc_diagnostics: the diagnostics do not match the stored points.")

  output <- data.frame(Chain = rep(seq_along(lengths), lengths),
                       Iteration = sequence(lengths))
  for (k in seq_along(names))
    output[[names[k]]] <- values[, k]
  # (in a scan, one set of columns per move: moveK.divergent, …)
  for (name in grep("(^|\\.)divergent$", names(output), value = TRUE))
    output[[name]] <- output[[name]] != 0
  output
}

#' Posterior expectations from a run with factor estimates that may be negative.
#'
#' A run whose factor estimates may be negative (normaliser = estimate(...,
#' debias = roulette(...), block_poisson(...) or taylor(...))) targets the
#' absolute value of the estimated posterior and records each point's sign.
#' The expectation of a function of the parameters under the posterior is
#' then the mean of the function weighted by the signs (Lyne et al. 2015):
#' sum(s_i f_i)/sum(s_i) for MCMC, and sum(W_i s_i f_i)/sum(W_i s_i) for
#' importance sampling and SMC, with W_i the points' normalised weights.
#'
#' @param output The output of load_mcmc_output, or of load_smc_output, with its Sign column.
#' @param parameter The name of a parameter.
#' @param dimension Its element (default 1).
#' @param f A function of the parameter's values (default: the identity).
#' @param burn_in MCMC: the number of first iterations of each chain to leave out.
#' @param target IS and SMC: the target whose points are used (default: the last).
#' @return A list: the estimate, and the (weighted) mean sign (near 0 makes the estimate unreliable).
#' @export
signed_mean <- function(output, parameter, dimension = 1, f = identity, burn_in = 0, target = NULL)
{
  if (!("Sign" %in% names(output)))
    stop("signed_mean: the output has no Sign column (no factor estimate of the run could be negative).")
  if ("LogWeight" %in% names(output) && "Target" %in% names(output))
  {
    # importance sampling and SMC: the points of one target, weighted
    if (is.null(target))
      target = max(output$Target)
    rows = output[(output$ParameterName == parameter) & (output$Dimension == dimension) & (output$Target == target), , drop = FALSE]
    if (nrow(rows) == 0)
      stop("signed_mean: no values of ", parameter, "[", dimension, "] at target ", target, ".")
    w = exp(rows$LogWeight - max(rows$LogWeight))
    w = w / sum(w)
    s = rows$Sign
    if (sum(w * s) == 0)
      stop("signed_mean: the weighted signs sum to zero, so the estimate is not defined.")
    return(list(estimate = sum(w * s * f(rows$Value)) / sum(w * s), mean_sign = sum(w * s)))
  }
  rows = output[(output$ParameterName == parameter) & (output$Dimension == dimension) & (output$Iteration > burn_in), , drop = FALSE]
  if (nrow(rows) == 0)
    stop("signed_mean: no values of ", parameter, "[", dimension, "] after the burn-in.")
  s = rows$Sign
  if (sum(s) == 0)
    stop("signed_mean: the signs sum to zero, so the estimate is not defined.")
  list(estimate = sum(s * f(rows$Value)) / sum(s), mean_sign = mean(s))
}
