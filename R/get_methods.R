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

      if ( !is.null(names(parameter_list)) && (h %in% names(parameter_list)) )
      {
        # a run-time constant given by name (dev/plan.md 3.14, L-6)
        parameter_arguments = c(parameter_arguments,h)
        do.call("<-",list(h, parameter_list[[h]]))
      }
      else if ( (nchar(h)>1) && (substr(h,1,1)=="p") && (!grepl("\\D",substr(h,2,nchar(h)))) )
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

# "name = value" at the top level of an argument: list(name, value text);
# otherwise list("", text)
split_keyword_argument = function(text)
{
  text = trimws(text)
  m = regmatches(text, regexec("^([A-Za-z_.][A-Za-z0-9_.]*)\\s*=([^=].*)$", text))[[1]]
  if (length(m) == 3)
    return(list(trimws(m[2]), trimws(m[3])))
  list("", text)
}

# The SMC sequences (/***smc_sequence,ilike::annealing(power,c(0,0.5,1))***/),
# or with method_name = "enk_sequence" the EnK sequences (EnKI),
# as one built-in call per sequence for the C++ registry (dev/plan.md 3.14,
# I2b): list(sequences = list(list(builtin, arguments), ...)). The variable is
# passed as text; the other arguments (the schedule, the factors) are R
# expressions, evaluated here with the run-time constants (by name, or as
# p1, p2, ...). The registry checks the names and numbers of arguments.
get_smc_sequences = function(model,model_parameters,method_name="smc_sequence")
{
  # no method blocks: no sequence (the caller uses its default)
  if (!("method" %in% names(model)))
    return(NULL)

  sequences = list()
  for (method in model[["method"]])
  {
    if (!(method_name %in% names(method)))
      next
    block = method[[method_name]]
    if (!all(c("type","variables","parameters") %in% names(block)))
      next

    # the variable, then the arguments: text from the old format (R
    # expressions, evaluated here), or values from the ilike language
    items = c(list(block[["variables"]]), block[["parameters"]])
    item_names = c("", if (is.null(names(block[["parameters"]]))) rep("", length(block[["parameters"]])) else names(block[["parameters"]]))
    arguments = list()
    argument_names = character(0)
    for (k in seq_along(items))
    {
      item = items[[k]]
      if (!is.character(item) || length(item) != 1)
      {
        arguments[[length(arguments)+1]] = item
        argument_names = c(argument_names, item_names[k])
        next
      }
      if (nchar(trimws(item)) == 0)
        next
      kv = if (item_names[k] != "") list(item_names[k], item) else split_keyword_argument(item)
      name = kv[[1]]
      if ( (name == "variable") || ((name == "") && (k == 1)) )
      {
        value = kv[[2]]
      }
      else if (item_names[k] != "")
      {
        # a name given in the ilike language: a constant
        value = if (!is.null(model_parameters[[kv[[2]]]])) model_parameters[[kv[[2]]]] else kv[[2]]
      }
      else
      {
        f_info = method_parse(kv[[2]],model_parameters)
        if (length(f_info[[2]])>0)
          stop(paste(method_name,' ilike::',block[["type"]],': ',kv[[2]],' uses ',paste(f_info[[2]],collapse=", "),
                     ', which is not a run-time constant (in model_parameter_list).',sep=""))
        value = f_info[[1]]()
      }
      arguments[[length(arguments)+1]] = value
      argument_names = c(argument_names, name)
    }
    names(arguments) = argument_names
    sequences[[length(sequences)+1]] = list(builtin = block[["type"]], arguments = arguments)
  }

  if (length(sequences) == 0)
    return(NULL)
  list(sequences = sequences)
}

# With no sequence given: anneal from the prior (power 0) to the posterior
# (power 1) in one step.
default_annealing_sequence = function()
{
  list(sequences = list(list(builtin = "annealing", arguments = list("power", c(0,1)))))
}

# The EnK sequences (/***enk_sequence,ilike::annealing(power,c(0,0.5,1))***/),
# read as the SMC sequences; the C++ registry allows only ilike::annealing.
get_enk_sequences = function(model,model_parameters)
{
  get_smc_sequences(model,model_parameters,"enk_sequence")
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

# TRUE if the recipe, or a recipe nested in it (an algorithmic likelihood's
# model), contains functions that call R (R expressions; R, Python or Julia
# functions).
recipe_uses_r_functions = function(recipe)
{
  if (!is.list(recipe))
    return(FALSE)
  if (isTRUE(recipe[["uses_r_functions"]]))
    return(TRUE)
  any(vapply(recipe, recipe_uses_r_functions, logical(1)))
}

# For parallel = TRUE in IS, SMC_with_MCMC and MCMC (dev/plan.md 2.6a, G4):
# refuse models that call R (R is single-threaded), and pass on the number of
# threads RcppParallel is set to use (RcppParallel::setThreadOptions).
prepare_parallel = function(recipe, parallel)
{
  # a new run: its start-of-run notes are printed once each (dev/plan.md 3.11)
  reset_informs()
  # the number of threads is passed even when parallel = FALSE: an inner
  # particle filter may be asked to run in parallel (ilike::pf(..., TRUE))
  threads = suppressWarnings(as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")))
  if (is.na(threads) || (threads < 1))
    threads = RcppParallel::defaultNumThreads()
  recipe$number_of_threads = threads
  if (isTRUE(parallel) && recipe_uses_r_functions(recipe))
    stop("parallel = TRUE needs a model written in C++ (or with ilike:: built-in functions): this model has functions written in R, Python or Julia, or R expressions, which cannot be called from parallel threads.")
  recipe
}
