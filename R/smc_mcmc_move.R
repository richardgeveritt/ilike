#' SMC with an MCMC move
#'
#' @param model A file containing the model, or a pre-compiled model list.
#' @param results_directory The name of the directory to which results will be written.
#' @param number_of_particles The number of particles.
#' @param number_of_mcmc_iterations (optional) The number of MCMC iterations.
#' @param mcmc_termination_method (optional) The method used to terminate the MCMC runs.
#' @param adaptive_resampling_ess (optional) Resample each time the ESS drops below this value.
#' @param adaptive_resampling_method (optional) Specify more generally a method to decide when to reasample.
#' @param adaptive_target_method (optional) Specify a method to decide how to adapt the sequence of targets.
#' @param smc_termination_method (optional) Method that determines when the SMC terminates (prior to finishing the sequence of targets).
#' @param mcmc_at_last_step (optional) Do we run the MCMC at the final step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param seed (optional) The seed for the random number generator.
#' @param parallel_flag (optional) Set to true to perform the importance sampling in parallel, false for serial.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Nothing: output can be found in the output_directory.
#' @export
smc_mcmc_move = function(model,
                         results_directory,
                         number_of_particles,
                         number_of_mcmc_iterations = 1,
                         mcmc_termination_method = NULL,
                         adaptive_resampling_ess = NULL,
                         adaptive_resampling_method = NULL,
                         adaptive_target_method = NULL,
                         smc_termination_method = NULL,
                         mcmc_at_last_step = FALSE,
                         smc_iterations_to_store = 2,
                         write_to_file_at_each_iteration = TRUE,
                         model_parameter_list = list(),
                         algorithm_parameter_list = list(),
                         seed = NULL,
                         parallel_flag = FALSE,
                         grain_size = 100000)
{
  if ((is.character(model)) && (length(model) == 1))
    model = parse_ilike_model(model,model_parameter_list)

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  # Sort MCMC termination method.
  if (!is.null(mcmc_termination_method))
  {
    print("MCMC termination method provided: using this method instead of number_of_mcmc_iterations.")
    if ("values" %in% names(mcmc_termination_method))
    {
      mcmc_termination_method[["values"]] = paste(mcmc_termination_method[["values"]])
    }
  }
  else
  {
    mcmc_termination_method = list(method="iterations",values=list(as.character(number_of_mcmc_iterations)))
  }

  # Adaptive resampling method.
  if (is.null(adaptive_resampling_ess) && is.null(adaptive_resampling_method) )
  {
    print("No method set for adaptive resampling: defaulting to resampling at every iteration.")
    adaptive_resampling_ess = number_of_particles
    adaptive_resampling_method = list(criterion='ess',values=paste(adaptive_resampling_ess))
  }
  else if (!is.null(adaptive_resampling_ess) && is.null(adaptive_resampling_method) )
  {
    print("No method set for adaptive resampling: using adaptive_resampling_ess.")
    adaptive_resampling_method = list(criterion='ess',values=paste(adaptive_resampling_ess))
  }
  else if (is.null(adaptive_resampling_ess) && !is.null(adaptive_resampling_method) )
  {
    if ("values" %in% names(adaptive_resampling_method))
    {
      adaptive_resampling_method[["values"]] = paste(adaptive_resampling_method[["values"]])
    }
  }
  else
  {
    print("Method set for adaptive resampling: ignoring adaptive_resampling_ess.")
    if ("values" %in% names(adaptive_resampling_method))
    {
      adaptive_resampling_method[["values"]] = paste(adaptive_resampling_method[["values"]])
    }
  }

  # Adaptive target.
  if (is.null(adaptive_target_method))
  {
    print("No method set for adaptive sequence of target distributions in SMC: defaulting to moving directly from one distribution to another.")
    adaptive_target_method = list(criterion='positive')
  }
  else
  {
    if ("values" %in% names(adaptive_target_method))
    {
      adaptive_target_method[["values"]] = paste(adaptive_target_method[["values"]])
    }
  }

  # SMC sequencing.
  if (is.null(smc_sequencer_method))
  {
    print("No SMC sequence set. Defaulting to importance sampling with prior as proposal.")
    smc_sequencer_method = list(types=c("annealing"),variables=c("power"),sequences=paste(c(0,1)))
  }
  else
  {
    if ("sequences" %in% names(smc_sequencer_method))
    {
      if (is.vector(smc_sequencer_method[["sequences"]]))
      {
        smc_sequencer_method[["sequences"]] = paste(smc_sequencer_method[["sequences"]])
      }
      else if (is.list(smc_sequencer_method[["sequences"]]))
      {
        for (i in 1:length(smc_sequencer_method[["sequences"]]))
        {
          if (is.vector(smc_sequencer_method[["sequences"]][[i]]))
          {
            smc_sequencer_method[["sequences"]][[i]] = paste(smc_sequencer_method[["sequences"]][[i]])
          }
        }
      }
      else
      {
        stop('smc_sequencer_method[["sequences"]] must be a numeric vector or a list of numeric vectors.')
      }
    }
  }

  # SMC termination.
  if (is.null(smc_termination_method))
  {
    smc_termination_method = list()
  }
  else
  {
    if ("values" %in% names(smc_termination_method))
    {
      smc_termination_method[["values"]] = paste(smc_termination_method[["values"]])
    }
  }

  do_smc_mcmc_move(model,
                   model_parameter_list,
                   algorithm_parameter_list,
                   number_of_particles,
                   mcmc_termination_method,
                   adaptive_resampling_method,
                   smc_sequencer_method,
                   smc_termination_method,
                   smc_iterations_to_store,
                   write_to_file_at_each_iteration,
                   parallel_flag,
                   grain_size,
                   results_directory,
                   seed)
}
