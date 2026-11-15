method_parse <- function(input,
                         parameter_list = list())
{
  # The variables used in the expression (function names and named-argument
  # names are not variables).
  parsed = tryCatch(parse(text=input, keep.source=FALSE),
                    error = function(e) stop(paste("Could not parse the R expression ", input, ": ", conditionMessage(e), sep="")))
  required_args <- unique(all.vars(parsed))

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
          stop(paste("In smc_sequence, parameter ",as.numeric(substr(h,2,nchar(h)))," not found in parameter_list (parameter_list has length ",length(parameter_list),").",sep=""))
        }
      }
    }

    # If p1, etc are in the list, remove from arg list, set p1=, etc. Then should be used by function.
    required_args = setdiff(required_args,parameter_arguments)
  }

  result_function = eval(parse(text = paste('result_function <- function(', paste(required_args,collapse=","), ') { return(' , input , ')}', sep='')))
  #result_function_text = paste(' <- function(', paste(required_args,collapse=","), ') { return(' , input , ')}', sep='')
  return(list(result_function,required_args))
}

get_method = function(model,method_name)
{
  output_method = NULL

  if ("method" %in% names(model))
  {
    methods = model[["method"]]
    for (i in 1:length(methods))
    {
      if ( (method_name=="mcmc_weights") || (method_name=="enk_likelihood_index") )
      {
        if (method_name %in% names(methods[[i]]))
        {
          output_method = list(func=methods[[i]][[method_name]])
        }
      }
      else
      {
        if (method_name %in% names(methods[[i]]) && ("type" %in% names(methods[[i]][[method_name]])) && ("parameters" %in% names(methods[[i]][[method_name]])) )
        {
          output_method = list(method=methods[[i]][[method_name]][["type"]],values=methods[[i]][[method_name]][["parameters"]])
        }
      }
    }
  }
  else
  {
    return(output_method)
  }

  return(output_method)
}

get_smc_sequences = function(model,model_parameters)
{
  method_name = "smc_sequence"

  output_method = NULL

  counter = 0

  if ("method" %in% names(model))
  {
    methods = model[["method"]]
    for (i in 1:length(methods))
    {
      if (method_name %in% names(methods[[i]]) && ("type" %in% names(methods[[i]][[method_name]])) && ("variables" %in% names(methods[[i]][[method_name]])) && ("parameters" %in% names(methods[[i]][[method_name]])) )
      {
        if (methods[[i]][[method_name]][["type"]]=="parameter")
        {
          sequence = methods[[i]][[method_name]][["parameters"]][[1]]

          f_info = method_parse(sequence,model_parameters)

          if (length(f_info[[2]])>0)
          {
            stop('SMC sequence must not have any variables (except for parameters "p1", "p2", etc.)')
          }

          # Optional third argument: the factors (1-based) that depend on the
          # sequence variable. If it is not given, every factor is
          # re-evaluated when the target changes.
          if (length(methods[[i]][[method_name]][["parameters"]])>2)
          {
            stop('smc_sequence ilike::parameter: too many arguments (expected the variable name, the sequence of values and optionally the factors affected).')
          }
          factors_affected = NULL
          if (length(methods[[i]][[method_name]][["parameters"]])==2)
          {
            factors_affected = method_parse(methods[[i]][[method_name]][["parameters"]][[2]],model_parameters)[[1]]()
          }

          if (counter==0)
          {
            schedule_list = list()
            schedule_list[[1]] = f_info[[1]]()
            output_method = list(types=methods[[i]][[method_name]][["type"]],variables=methods[[i]][[method_name]][["variables"]],schedules=schedule_list)
          }
          else
          {
            output_method[["types"]] = c(output_method[["types"]],methods[[i]][[method_name]][["type"]])
            output_method[["variables"]] = c(output_method[["variables"]],methods[[i]][[method_name]][["variables"]])
            output_method[["schedules"]] = append(output_method[["schedules"]],f_info[[1]]())
          }
          if (!is.null(factors_affected))
          {
            output_method[["factors_affected"]] = append(output_method[["factors_affected"]],factors_affected)
          }
        }
        else
        {
          # Annealing/tempering: (variable, schedule) or (variable, schedule,
          # factors). The optional factors (1-based) are the ones raised to the
          # power; the others are not annealed.
          sequence_parameters = methods[[i]][[method_name]][["parameters"]]
          if (length(sequence_parameters)>2)
          {
            stop(paste('smc_sequence ilike::',methods[[i]][[method_name]][["type"]],': too many arguments (expected the variable name, the schedule and optionally the factors to anneal).',sep=""))
          }
          sequence = sequence_parameters[[1]]
          factors_affected = NULL
          if (length(sequence_parameters)==2)
          {
            factors_affected = method_parse(sequence_parameters[[2]],model_parameters)[[1]]()
          }

          f_info = method_parse(sequence,model_parameters)

          if (length(f_info[[2]])>0)
          {
            stop('SMC sequence must not have any variables (except for parameters "p1", "p2", etc.)')
          }

          if (counter==0)
          {
            schedule_list = list()
            schedule_list[[1]] = f_info[[1]]()
            output_method = list(types=methods[[i]][[method_name]][["type"]],variables=methods[[i]][[method_name]][["variables"]],schedules=schedule_list)
          }
          else
          {
            output_method[["types"]] = c(output_method[["types"]],methods[[i]][[method_name]][["type"]])
            output_method[["variables"]] = c(output_method[["variables"]],methods[[i]][[method_name]][["variables"]])
            output_method[["schedules"]] = append(output_method[["schedules"]],f_info[[1]]())
          }
          if (!is.null(factors_affected))
          {
            output_method[["factors_affected"]] = append(output_method[["factors_affected"]],factors_affected)
          }
        }

        counter = counter + 1

      }
    }
  }
  else
  {
    stop(paste("Model file needs to specify a valid method for ",method_name,".",sep=""))
  }

  return(output_method)
}

get_enk_sequences = function(model,model_parameters)
{
  method_name = "enk_sequence"

  output_method = NULL

  counter = 0

  if ("method" %in% names(model))
  {
    methods = model[["method"]]
    for (i in 1:length(methods))
    {
      if (method_name %in% names(methods[[i]]) && ("type" %in% names(methods[[i]][[method_name]])) && ("variables" %in% names(methods[[i]][[method_name]])) && ("parameters" %in% names(methods[[i]][[method_name]])) )
      {
        sequence = methods[[i]][[method_name]][["parameters"]]

        f_info = method_parse(sequence,model_parameters)

        if (length(f_info[[2]])>0)
        {
          stop('SMC sequence must not have any variables (except for parameters "p1", "p2", etc.')
        }

        if (counter==0)
        {
          schedule_list = list()
          schedule_list[[1]] = f_info[[1]]()
          output_method = list(types=methods[[i]][[method_name]][["type"]],variables=methods[[i]][[method_name]][["variables"]],schedules=schedule_list)
        }
        else
        {
          output_method[["types"]] = c(output_method[["types"]],methods[[i]][[method_name]][["type"]])
          output_method[["variables"]] = c(output_method[["variables"]],methods[[i]][[method_name]][["variables"]])
          output_method[["schedules"]] = append(output_method[["schedules"]],f_info[[1]]())
        }

        counter = counter + 1

      }
    }
  }
  else
  {
    stop(paste("Model file needs to specify a valid method for ",method_name,".",sep=""))
  }

  return(output_method)
}

# The parallel workers (RcppParallelSMCWorker and the ensemble Kalman
# equivalents) are disconnected in the C++ code: with parallel = TRUE some
# algorithms stop with "Parallel worker not set up" and others leave the worker
# NULL and crash. Refuse the option here until the workers are restored
# (dev/plan.md section 3.1).
stop_if_parallel = function(parallel)
{
  if (isTRUE(parallel))
    stop("parallel = TRUE is not available yet: run with parallel = FALSE.")
}
