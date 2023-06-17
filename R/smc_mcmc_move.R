#' SMC with an MCMC move
#'
#' @param model A file containing the model, or a pre-compiled model list.
#' @param number_of_particles The number of particles.
#' @param parallel_flag (optional) Set to true to perform the importance sampling in parallel, false for serial.
#' @param results_directory (optional) The name of the directory to which results will be written.
#' @param smc_iterations_to_store (optional) The number of iterations of SMC output stored in memory as the algorithm is running (cannot be fewer than 2).
#' @param write_to_file_at_each_iteration (optional) Do we write the algorithm output to file at each SMC step (TRUE/FALSE)?
#' @param model_parameter_list (optional) A list containing parameters for the model.
#' @param algorithm_parameter_list (optional) A list containing named parameters for the algorithm.
#' @param seed (optional) The seed for the random number generator.
#' @param grain_size (optional) Sets a minimum chunk size for parallelisation (see https://oneapi-src.github.io/oneTBB/main/tbb_userguide/Controlling_Chunking_os.html).
#' @return Estimate of the marginal likelihood.
#' @export
smc_mcmc_move = function(model,
                         number_of_particles,
                         parallel_flag = FALSE,
                         results_directory = getwd(),
                         smc_iterations_to_store = 2,
                         write_to_file_at_each_iteration = TRUE,
                         model_parameter_list = list(),
                         algorithm_parameter_list = list(),
                         seed = NULL,
                         grain_size = 100000)
{
  if ((is.character(model)) && (length(model) == 1))
    model = parse_ilike_model(model,model_parameter_list)

  if (is.null(seed))
  {
    seed = ilike_rdtsc()
  }

  # Sort MCMC termination method.
  mcmc_termination_method = get_method(model,"mcmc_termination")

  if (is.null(mcmc_termination_method))
  {
    stop("No MCMC termination method provided: using default of 10^5 iterations.")
    mcmc_termination_method = list(method="iterations",values=list(as.character(100000)))
  }

  # Adaptive resampling method.
  adaptive_resampling_method = get_method(model,"adaptive_resampling")

  if (is.null(adaptive_resampling_method))
  {
    print("No method set for adaptive resampling: using adaptive_resampling_ess.")
    adaptive_resampling_method = list(criterion='ess',values=as.list(paste(adaptive_resampling_ess)))
  }

  # Adaptive target method.
  adaptive_target_method = get_method(model,"adaptive_target")

  if (is.null(adaptive_target_method))
  {
    print("No method set for adaptive sequence of target distributions in SMC: defaulting to moving directly from one distribution to another.")
    adaptive_target_method = list(criterion='positive')
  }

  # SMC sequencing.
  smc_sequencer_method = get_smc_sequences(model,model_parameter_list)

  if (is.null(smc_sequencer_method))
  {
    print("No SMC sequence set. Defaulting to importance sampling with prior as proposal.")
    smc_sequencer_method = list(types=c("annealing"),variables=c("power"),schedules=c(0,1))
  }

  # SMC termination method.
  smc_termination_method = get_method(model,"smc_termination")

  if (is.null(smc_termination_method))
  {
    smc_termination_method = list()
  }

  return(do_smc_mcmc_move(model,
                          model_parameter_list,
                          algorithm_parameter_list,
                          number_of_particles,
                          mcmc_termination_method,
                          adaptive_resampling_method,
                          smc_sequencer_method,
                          adaptive_target_method,
                          smc_termination_method,
                          smc_iterations_to_store,
                          write_to_file_at_each_iteration,
                          parallel_flag,
                          grain_size,
                          results_directory,
                          seed))
}