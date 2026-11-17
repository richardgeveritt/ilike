# Converting recipes from the old .ilike format to the ilike language
# (dev/ilike-language.md section 5; dev/plan.md 3.14, I3c). The old parser
# reads the files in a record-only mode, so that the blocks are grouped and
# numbered exactly as compile() groups and numbers them; the groups are then
# written in the language. The C++ code is copied unchanged.

#' Convert a recipe from the old .ilike format to the ilike language
#'
#' Rewrites the files of a recipe written in the old format (with /***...***/
#' blocks) as one file in the ilike language (see the package's documentation of
#' the language). The files are read as compile() reads them: blocks are grouped
#' into factors, proposals and settings by the same rules, factors are named
#' (a built-in prior by its variable, others f1, f2, ... by their number) and
#' the old run-time constants p1, p2, ... are declared as constants of those
#' names (an unnamed model_parameter_list is still matched to them in order).
#' All moves and settings go into one algorithm section called main. The inner
#' file of an algorithmic likelihood (ilike::pf("inner.ilike", ...)) is
#' converted too, into the folder of the output, and imported.
#'
#' R expressions become calls in the language (r::f(...), or py::/jl:: for
#' functions in the recipe's Python or Julia file, which is copied next to the
#' output and imported); an expression the language cannot read directly
#' becomes an R function in an r section. Recipes with potential functions,
#' proposal transforms, EnK transforms or enk_likelihood_index blocks cannot be
#' converted yet.
#'
#' @param filenames The file(s) of the recipe, as given to compile().
#' @param output The file to write.
#' @param overwrite (optional) If TRUE, an existing output file is replaced.
#' @return Invisibly, the converted text.
#' @export
convert_ilike_file <- function(filenames, output, overwrite = FALSE)
{
  for (f in filenames)
    if (!file.exists(f))
      stop(paste("File", f, "does not exist."))
  if (file.exists(output) && !overwrite)
    stop(paste("The output file ", output, " exists: use overwrite = TRUE to replace it.", sep = ""))
  if (any(normalizePath(filenames) == normalizePath(output, mustWork = FALSE)))
    stop("The output file must not be one of the input files.")
  for (f in filenames)
    if (is_ilike_language_text(paste(readLines(f, warn = FALSE), collapse = "\n")))
      stop(paste("The file ", f, " is already in the ilike language.", sep = ""))

  # the old parser, recording the blocks instead of compiling them
  .ilike_env$recorded_blocks = list()
  on.exit(.ilike_env$recorded_blocks <- NULL, add = TRUE)
  records = compile(filenames)
  .ilike_env$recorded_blocks = NULL

  # helper files with the recipe files' base names (.R, .py, .jl), which the
  # old format loads: copied next to the output and imported
  .ilike_env$convert_state = list(r_code = character(0), count = 0, imports = character(0),
                                  py_functions = character(0), jl_functions = character(0))
  on.exit(.ilike_env$convert_state <- NULL, add = TRUE)
  for (f in filenames)
  {
    stem = tools::file_path_sans_ext(f)
    for (extension in c("R", "py", "jl"))
    {
      helper = paste(stem, ".", extension, sep = "")
      if (!file.exists(helper))
        next
      copy = file.path(dirname(output), basename(helper))
      if (normalizePath(copy, mustWork = FALSE) != normalizePath(helper))
        file.copy(helper, copy, overwrite = TRUE)
      alias = c(R = "r", py = "py", jl = "jl")[[extension]]
      .ilike_env$convert_state$imports = c(.ilike_env$convert_state$imports,
                                           paste("import \"", basename(helper), "\" as ", alias, sep = ""))
      lines = readLines(helper, warn = FALSE)
      if (extension == "py")
        .ilike_env$convert_state$py_functions = c(.ilike_env$convert_state$py_functions,
                                                  sub("^def\\s+(\\w+).*$", "\\1", grep("^def\\s+\\w+", lines, value = TRUE)))
      if (extension == "jl")
        .ilike_env$convert_state$jl_functions = c(.ilike_env$convert_state$jl_functions,
                                                  sub("^function\\s+(\\w+).*$", "\\1", grep("^function\\s+\\w+", lines, value = TRUE)))
    }
  }

  text = convert_records(records, filenames, output)
  writeLines(text, output)
  invisible(text)
}

# --- reading the old blocks ---------------------------------------------------

convert_stop <- function(record, message)
{
  stop(paste("convert_ilike_file: line ", record$line, " (block ", record$name, "): ", message, sep = ""), call. = FALSE)
}

# "ilike::name(a, b)" as list(name, arguments as text)
convert_builtin_call <- function(record)
{
  call = gsub("\\s+", "", record$block_function)
  if (!grepl("^ilike::", call))
    convert_stop(record, "this block is an R expression where the language needs a built-in (ilike::...).")
  call = sub("^ilike::", "", call)
  name = sub("\\(.*$", "", call)
  inner = substr(call, nchar(name) + 2, nchar(call) - 1)
  arguments = if (nchar(inner) == 0) character(0) else split_string_at_comma_ignoring_parentheses(inner)
  list(name = name, arguments = arguments[nchar(trimws(arguments)) > 0])
}

# an argument of a built-in, from the old format's text to the language's
convert_argument <- function(text, record)
{
  t = trimws(text)
  m = regmatches(t, regexec("^([[:alpha:]_][[:alnum:]_]*)\\s*=([^=].*)$", t))[[1]]
  if (length(m) == 3)
  {
    name = if (m[2] == "length.out") "length" else m[2]
    return(paste(name, "=", convert_argument(m[3], record)))
  }
  if (grepl("^c\\(.*\\)$", t))
  {
    items = split_string_at_comma_ignoring_parentheses(substr(t, 3, nchar(t) - 1))
    if (length(items) == 1)
      return(convert_argument(items, record)) # c(x) is x
    return(paste("[", paste(sapply(items, convert_argument, record = record), collapse = ", "), "]", sep = ""))
  }
  if (grepl("^(seq|rep|diag)\\(.*\\)$", t))
  {
    name = sub("\\(.*$", "", t)
    items = split_string_at_comma_ignoring_parentheses(substr(t, nchar(name) + 2, nchar(t) - 1))
    return(paste(name, "(", paste(sapply(items, convert_argument, record = record), collapse = ", "), ")", sep = ""))
  }
  if (grepl("^'.*'$", t))
    return(paste("\"", substr(t, 2, nchar(t) - 1), "\"", sep = ""))
  if (grepl("^\".*\"$", t))
    return(t)
  if (grepl("^[-+]?([0-9.]+([eE][-+]?[0-9]+)?|Inf|TRUE|FALSE)$", t) || grepl("^[[:alpha:]_][[:alnum:]_]*$", t))
    return(t)
  if (grepl("^[-+*/^ ().0-9eE]+$", t))
    return(t)
  convert_stop(record, paste("the argument ", t, " is an R expression the language does not have.", sep = ""))
}

convert_call_text <- function(name, arguments, record)
{
  paste(name, "(", paste(sapply(arguments, convert_argument, record = record), collapse = ", "), ")", sep = "")
}

# the name of a block's C++ function (as the old parser finds it)
convert_function <- function(record)
{
  if (!isTRUE(record$custom))
  {
    foreign = convert_foreign(record)
    if (!is.null(foreign$output) && record$name != "data")
      convert_stop(record, "unexpected output variable.")
    return(foreign$call)
  }
  paste("cpp::", trimws(function_name_for_cpp(record$code)), sep = "")
}

# An R expression of the old format (which may call R, Python or Julia
# functions) as a call in the language: r::f(...) (py::, jl:: when f is in
# the recipe's Python or Julia file) if it is one call the language can
# read; otherwise an R function in an r section, called with the objects the
# expression uses. "x = expression" (simulating, or the data) gives a named
# list (simulating) or the data variable x.
convert_foreign <- function(record)
{
  text = trimws(record$block_function)
  # the documented parameters.tau form, as parameters$tau
  text = gsub("\\b(parameters|data|proposed_parameters|proposal_parameters)\\.([^[:space:](),=+*/^<>!&|-]+)", "\\1$\\2", text, perl = TRUE)
  output = NULL
  position = find_top_level_assignment(text)
  if (position > 0)
  {
    output = trimws(substr(text, 1, position - 1))
    text = trimws(substr(text, position + 1, nchar(text)))
    if (record$name != "data")
    {
      text = paste("list(", output, " = ", text, ")", sep = "")
      output = NULL
    }
  }
  state = .ilike_env$convert_state
  function_name = sub("\\(.*$", "", text)
  readable = grepl("^[A-Za-z_][A-Za-z0-9_.]*\\(.*\\)$", text) &&
    tryCatch({ ilike_language_syntax_tree(paste("ilike 1.0\nmodel { factor x = r::", text, " }", sep = ""), "", FALSE); TRUE },
             error = function(e) FALSE) &&
    !grepl("\\.", function_name)
  # (only if its free names are the objects a role passes or constants p1, ...:
  # in the language a bare name is a parameter, but in R it may be anything,
  # e.g. pi)
  objects = c("variable", "variable1", "variable2", "proposed_parameters", "parameters", "proposal_parameters", "data")
  if (readable)
  {
    free = all.vars(parse(text = text))
    members = regmatches(text, gregexpr("\\$[^[:space:](),=+*/^<>!&|$-]+", text))[[1]]
    free = setdiff(free, c(objects, sub("^\\$", "", members)))
    readable = all(grepl("^p[0-9]+$", free))
  }
  if (readable)
  {
    name_space = if (function_name %in% state$py_functions) "py" else if (function_name %in% state$jl_functions) "jl" else "r"
    return(list(output = output, call = paste(name_space, "::", text, sep = "")))
  }
  # an R function, called with the objects the expression uses
  used = objects[sapply(objects, function(o) grepl(paste("\\b", o, "\\b", sep = ""), text, perl = TRUE))]
  state$count = state$count + 1
  name = paste("ilike_converted_", state$count, sep = "")
  state$r_code = c(state$r_code, paste("# from ", record$name, " (line ", record$line, ")", sep = ""),
                   paste(name, " = function(", paste(used, collapse = ", "), ") {", sep = ""),
                   paste("  ", text, sep = ""), "}")
  .ilike_env$convert_state = state
  list(output = output, call = paste("r::", name, "(", paste(used, collapse = ", "), ")", sep = ""))
}

# the string a small C++ function returns (return "x";)
convert_returned_string <- function(record)
{
  # /***data_variable,"x"***/: the name given as a string
  literal = trimws(record$block_function)
  if (!isTRUE(record$custom) && grepl("^([\"']).*\\1$", literal))
    return(substr(literal, 2, nchar(literal) - 1))
  m = regmatches(record$code, regexec("return\\s*\"([^\"]*)\"\\s*;", record$code))[[1]]
  if (length(m) != 2)
    convert_stop(record, "the function does not simply return a string (return \"x\";), so the variable's name cannot be read.")
  m[2]
}

# factor numbers from c(1, 2) (an R expression) or NumericVector::create(1, 2)
convert_factor_numbers <- function(record)
{
  source = if (isTRUE(record$custom)) record$code else gsub("\\s+", "", record$block_function)
  m = regmatches(source, regexec("(c|NumericVector::create)\\(([0-9, ]*)\\)", source))[[1]]
  if (length(m) != 3)
    convert_stop(record, "the factor numbers cannot be read (expected c(1, 2) or NumericVector::create(1, 2)).")
  as.numeric(trimws(strsplit(m[3], ",")[[1]]))
}

# --- writing the language ---------------------------------------------------------

convert_records <- function(records, filenames, output)
{
  code = character(0)
  groups = list()
  group_keys = character(0)
  for (r in records)
  {
    if (isTRUE(r$custom) && nchar(trimws(r$code)) > 0)
      code = c(code, sub("^\n", "", r$code))
    if (r$type == "pre" || r$number <= 0)
      next
    key = paste(r$type, r$number)
    if (!(key %in% group_keys))
    {
      group_keys = c(group_keys, key)
      groups[[key]] = list(type = r$type, number = r$number, entries = list())
    }
    if (!is.null(groups[[key]]$entries[[r$name]]))
      convert_stop(r, "a second block of this kind in the same group.")
    groups[[key]]$entries[[r$name]] = r
  }

  unsupported = c("potential_function", "enk_transform")
  for (g in groups)
    if (g$type %in% unsupported)
      convert_stop(g$entries[[1]], paste(g$type, "blocks are not in the ilike language yet."))

  # the factors' names, by number
  factor_names = character(0)
  for (g in groups)
    if (g$type == "factor")
    {
      name = paste("f", g$number, sep = "")
      if (!is.null(g$entries$prior))
      {
        variable = trimws(convert_builtin_call(g$entries$prior)$arguments[1])
        if (!(variable %in% factor_names))
          name = variable
      }
      factor_names[g$number] = name
    }
  name_of_factors = function(numbers, record)
  {
    if (any(numbers < 1 | numbers != round(numbers)))
      convert_stop(record, paste("the factor numbers ", paste(numbers, collapse = ", "), " are not factor numbers.", sep = ""))
    names = sapply(numbers, function(k) if (k <= length(factor_names) && !is.na(factor_names[k])) factor_names[k] else paste("f", k, sep = ""))
    if (any(numbers > length(factor_names) | is.na(factor_names[pmin(numbers, max(1, length(factor_names)))])))
      warning(paste("convert_ilike_file: line ", record$line, ": factors ", paste(numbers, collapse = ", "),
                    " are not in these files; they are written as ", paste(names, collapse = ", "),
                    ", which must be the factors' names in the model they are compiled with.", sep = ""), call. = FALSE)
    paste("[", paste(names, collapse = ", "), "]", sep = "")
  }

  data_lines = character(0)
  model_lines = character(0)
  algorithm_lines = character(0)
  import_lines = character(0)
  move_count = 0
  inner_count = 0

  for (g in groups)
  {
    e = g$entries
    if (g$type == "data")
    {
      if (length(data_lines) > 0)
        convert_stop(e[[1]], "a second data block (only the first is used).")
      if (!isTRUE(e$data$custom))
      {
        foreign = convert_foreign(e$data)
        data_lines = c(data_lines, paste("  ", if (is.null(foreign$output)) "data" else foreign$output, " = ", foreign$call, sep = ""))
      }
      else
        data_lines = c(data_lines, paste("  data =", convert_function(e$data)))
    }
    else if (g$type == "factor")
    {
      name = factor_names[g$number]
      line = convert_factor(e, name, filenames, output, function(path) {
        inner_count <<- inner_count + 1
        alias = paste("inner", inner_count, sep = "")
        import_lines <<- c(import_lines, paste("import \"", path, "\" as ", alias, sep = ""))
        alias
      })
      model_lines = c(model_lines, line)
    }
    else if (g$type == "transition_model")
      model_lines = c(model_lines, convert_transition(e, g$number))
    else if (g$type %in% c("m_proposal", "mh_proposal", "imh_proposal", "unadjusted_proposal"))
    {
      move_count = move_count + 1
      algorithm_lines = c(algorithm_lines, convert_move(e, g$type, paste("m", move_count, sep = ""), name_of_factors))
    }
    else if (g$type == "is_proposal")
    {
      if (!is.null(e$is_proposal))
      {
        b = convert_builtin_call(e$is_proposal)
        algorithm_lines = c(algorithm_lines, paste("  proposal =", convert_call_text(b$name, b$arguments, e$is_proposal)))
      }
      else
        algorithm_lines = c(algorithm_lines, paste("  proposal = density(simulate = ", convert_function(e$simulate_is_proposal),
                                                   ", log = ", convert_function(e$evaluate_log_is_proposal), ")", sep = ""))
    }
    else if (g$type == "transition_proposal")
    {
      if (!is.null(e$transition_proposal))
      {
        b = convert_builtin_call(e$transition_proposal)
        algorithm_lines = c(algorithm_lines, paste("  transition_proposal =", convert_call_text(b$name, b$arguments, e$transition_proposal)))
      }
      else
        algorithm_lines = c(algorithm_lines, paste("  transition_proposal = density(simulate = ", convert_function(e$simulate_transition_proposal),
                                                   ", log = ", convert_function(e$evaluate_log_transition_proposal), ")", sep = ""))
    }
    else if (g$type == "method")
    {
      out = convert_method(e, name_of_factors)
      if (out$where == "model") model_lines = c(model_lines, out$line) else algorithm_lines = c(algorithm_lines, out$line)
    }
    else
      convert_stop(e[[1]], paste("blocks of type", g$type, "are not in the ilike language yet."))
  }

  import_lines = c(.ilike_env$convert_state$imports, import_lines)
  body = c(import_lines, if (length(import_lines)) "" else NULL)
  all_text = paste(c(model_lines, algorithm_lines), collapse = "\n")
  constants = unique(regmatches(all_text, gregexpr("\\bp[0-9]+\\b", all_text))[[1]])
  constants = constants[order(as.numeric(sub("p", "", constants)))]
  if (length(constants) > 0)
    body = c(body, "constants {",
             paste("  constant ", constants, "   // the old format's run-time constant ", constants,
                   " (by name, or by position in an unnamed model_parameter_list)", sep = ""),
             "}", "")
  if (length(data_lines) > 0)
    body = c(body, "data {", data_lines, "}", "")
  if (length(model_lines) > 0)
    body = c(body, "model {", model_lines, "}", "")
  if (length(algorithm_lines) > 0)
    body = c(body, "algorithm main {", algorithm_lines, "}", "")
  if (length(.ilike_env$convert_state$r_code) > 0)
    body = c(body, "r {", .ilike_env$convert_state$r_code, "}", "")
  if (length(code) > 0)
  {
    # (runs of blank lines between the old blocks collapsed)
    code = strsplit(paste(code, collapse = "\n"), "\n")[[1]]
    blank = !nzchar(trimws(code))
    code = code[!(blank & c(FALSE, blank[-length(blank)]))]
    body = c(body, "cpp {", code, "}")
  }

  c("ilike 1.0",
    paste("// Converted by convert_ilike_file() from ", paste(basename(filenames), collapse = ", "),
          " (old format). The C++ code is unchanged.", sep = ""),
    "",
    body)
}

convert_factor <- function(e, name, filenames, output, import_inner)
{
  if (!is.null(e$prior))
  {
    b = convert_builtin_call(e$prior)
    variable = trimws(b$arguments[1])
    if (name == variable)
      return(paste("  ", variable, " ~ ", convert_call_text(b$name, b$arguments[-1], e$prior), sep = ""))
    return(paste("  prior ", name, " = ", convert_call_text(b$name, b$arguments, e$prior), sep = ""))
  }
  density = function(prefix_log, keys)
  {
    parts = character(0)
    for (k in names(keys))
      if (!is.null(e[[k]]))
        parts = c(parts, paste(keys[[k]], "=", convert_function(e[[k]])))
    if (length(parts) == 1 && grepl("^log =", parts))
      return(sub("^log = ", "", parts))
    paste("density(", paste(parts, collapse = ", "), ")", sep = "")
  }
  if (!is.null(e$evaluate_log_prior) || !is.null(e$simulate_prior))
    return(paste("  prior ", name, " = ",
                 density("prior", list(evaluate_log_prior = "log", simulate_prior = "simulate",
                                       evaluate_gradient_log_prior = "gradient",
                                       evaluate_second_gradient_log_prior = "second_gradient")), sep = ""))
  if (!is.null(e$evaluate_log_likelihood))
    return(paste("  factor ", name, " = ",
                 density("likelihood", list(evaluate_log_likelihood = "log",
                                            evaluate_gradient_log_likelihood = "gradient",
                                            evaluate_second_gradient_log_likelihood = "second_gradient")), sep = ""))
  if (!is.null(e$sbi_likelihood))
  {
    b = convert_builtin_call(e$sbi_likelihood)
    parts = c(convert_call_text(b$name, b$arguments, e$sbi_likelihood),
              paste("simulate =", convert_function(e$simulate_data_model)),
              paste("data =", gsub(";", ", ", convert_returned_string(e$data_variable))))
    if (grepl(",", parts[3]))
      parts[3] = sub("data = (.*)$", "data = [\\1]", parts[3])
    if (!is.null(e$summary_statistics))
      parts = c(parts, paste("summaries =", convert_function(e$summary_statistics)))
    return(paste("  factor ", name, " = sbi(", paste(parts, collapse = ", "), ")", sep = ""))
  }
  if (!is.null(e$linear_gaussian_data_model))
  {
    b = convert_builtin_call(e$linear_gaussian_data_model)
    a = sapply(b$arguments, convert_argument, record = e$linear_gaussian_data_model)
    return(paste("  factor ", name, " = linear_gaussian(data = ", a[1], ", state = ", a[2], ", matrix = ", a[3],
                 ", covariance = ", a[4], ")", sep = ""))
  }
  if (!is.null(e$linear_gaussian_data_matrix))
    return(paste("  factor ", name, " = linear_gaussian(data = ", convert_returned_string(e$linear_gaussian_data_variable),
                 ", state = ", convert_state_names(convert_returned_string(e$linear_gaussian_data_state_variable)),
                 ", matrix = ", convert_function(e$linear_gaussian_data_matrix),
                 ", covariance = ", convert_function(e$linear_gaussian_data_covariance), ")", sep = ""))
  if (!is.null(e$nonlinear_gaussian_data_function))
    return(paste("  factor ", name, " = nonlinear_gaussian(data = ", convert_returned_string(e$nonlinear_gaussian_data_variable),
                 ", function = ", convert_function(e$nonlinear_gaussian_data_function),
                 ", covariance = ", convert_function(e$nonlinear_gaussian_data_covariance), ")", sep = ""))
  if (!is.null(e$algorithmic_likelihood))
  {
    r = e$algorithmic_likelihood
    b = convert_builtin_call(r)
    # the inner file: converted into the output's folder, and imported
    inner = gsub("^[\"']|[\"']$", "", trimws(b$arguments[1]))
    candidates = c(inner, file.path(dirname(filenames), inner))
    found = candidates[file.exists(candidates)]
    if (length(found) == 0)
      convert_stop(r, paste("the inner file ", inner, " was not found.", sep = ""))
    inner_output = file.path(dirname(output), basename(found[1]))
    if (normalizePath(inner_output, mustWork = FALSE) == normalizePath(found[1]))
      convert_stop(r, "the inner file would be converted over itself: write the output in another folder.")
    if (!is_ilike_language_text(paste(readLines(found[1], warn = FALSE), collapse = "\n")))
      convert_ilike_file(found[1], inner_output, overwrite = TRUE)
    else
      file.copy(found[1], inner_output, overwrite = TRUE)
    alias = import_inner(basename(found[1]))
    rest = if (length(b$arguments) > 1) sapply(b$arguments[-1], convert_argument, record = r) else character(0)
    return(paste("  factor ", name, " = ", b$name, "(", paste(c(alias, rest), collapse = ", "), ")", sep = ""))
  }
  convert_stop(e[[1]], "this kind of factor is not in the ilike language yet.")
}

# "x;y" (several state variables) as [x, y]
convert_state_names <- function(names)
{
  parts = strsplit(names, ";")[[1]]
  if (length(parts) == 1) parts else paste("[", paste(parts, collapse = ", "), "]", sep = "")
}

convert_transition <- function(e, number)
{
  if (!is.null(e$transition_model))
  {
    b = convert_builtin_call(e$transition_model)
    return(paste("  transition ", trimws(b$arguments[1]), " = ", convert_call_text(b$name, b$arguments, e$transition_model), sep = ""))
  }
  if (!is.null(e$linear_gaussian_transition_model))
  {
    b = convert_builtin_call(e$linear_gaussian_transition_model)
    a = sapply(b$arguments, convert_argument, record = e$linear_gaussian_transition_model)
    return(paste("  transition ", a[1], " = linear_gaussian(matrix = ", a[2], ", covariance = ", a[3], ")", sep = ""))
  }
  if (!is.null(e$linear_gaussian_transition_matrix))
    return(paste("  transition ", convert_returned_string(e$linear_gaussian_transition_variable),
                 " = linear_gaussian(matrix = ", convert_function(e$linear_gaussian_transition_matrix),
                 ", covariance = ", convert_function(e$linear_gaussian_transition_covariance), ")", sep = ""))
  if (!is.null(e$nonlinear_gaussian_transition_function))
    return(paste("  transition ", convert_returned_string(e$nonlinear_gaussian_transition_variable),
                 " = nonlinear_gaussian(function = ", convert_function(e$nonlinear_gaussian_transition_function),
                 ", covariance = ", convert_function(e$nonlinear_gaussian_transition_covariance), ")", sep = ""))
  if (!is.null(e$simulate_transition_model))
  {
    parts = paste("simulate =", convert_function(e$simulate_transition_model))
    if (!is.null(e$evaluate_log_transition_model))
      parts = c(parts, paste("log =", convert_function(e$evaluate_log_transition_model)))
    return(paste("  transition t", number, " = density(", paste(parts, collapse = ", "), ")", sep = ""))
  }
  convert_stop(e[[1]], "this kind of transition model is not in the ilike language yet.")
}

convert_move <- function(e, type, move_name, name_of_factors)
{
  p = sub("_proposal$", "", type)
  kind = c(m = "metropolis", mh = "metropolis_hastings", imh = "independent_mh", unadjusted = "unadjusted")[[p]]
  for (k in names(e))
    if (grepl("transform", k))
      convert_stop(e[[k]], "proposal transforms are not in the ilike language yet.")
  variable = ""
  if (!is.null(e[[paste(p, "_proposal", sep = "")]]))
  {
    r = e[[paste(p, "_proposal", sep = "")]]
    b = convert_builtin_call(r)
    variable = trimws(b$arguments[1])
    parts = convert_call_text(b$name, b$arguments, r)
  }
  else
  {
    simulate = convert_function(e[[paste("simulate_", p, "_proposal", sep = "")]])
    log_record = e[[paste("evaluate_log_", p, "_proposal", sep = "")]]
    if (p == "imh")
      parts = paste("density(simulate = ", simulate, ", log = ", convert_function(log_record), ")", sep = "")
    else if (!is.null(log_record))
      parts = c(simulate, paste("log =", convert_function(log_record)))
    else
      parts = simulate
  }
  index = e[[paste(p, "_factor_index", sep = "")]]
  if (!is.null(index))
    parts = c(parts, paste("factors =", name_of_factors(convert_factor_numbers(index), index)))
  cut = e[[paste(p, "_cut_factor_index", sep = "")]]
  if (!is.null(cut))
    parts = c(parts, paste("cut_factors =", name_of_factors(convert_factor_numbers(cut), cut)))
  adaptor = e[[paste(p, "_adaptor", sep = "")]]
  adapt_function = e[[paste("adapt_", p, "_proposal", sep = "")]]
  if (!is.null(adapt_function))
  {
    f = convert_function(adapt_function)
    averaged = !is.null(adaptor) && convert_builtin_call(adaptor)$name == "averaged"
    parts = c(parts, paste("adapt = ", if (averaged) paste("averaged(", f, ")", sep = "") else f, sep = ""))
  }
  else if (!is.null(adaptor))
  {
    b = convert_builtin_call(adaptor)
    arguments = b$arguments
    # the adaptor's variable is the proposal's by default
    if (b$name != "acceptance_rate" && length(arguments) > 0 && trimws(arguments[1]) == variable)
      arguments = arguments[-1]
    parts = c(parts, paste("adapt =", convert_call_text(b$name, arguments, adaptor)))
  }
  if (!is.null(e$smc_adaptor))
  {
    b = convert_builtin_call(e$smc_adaptor)
    arguments = b$arguments
    if (length(arguments) > 0 && trimws(arguments[1]) == variable)
      arguments = arguments[-1]
    parts = c(parts, paste("smc_adapt =", convert_call_text(b$name, arguments, e$smc_adaptor)))
  }
  paste("  move ", move_name, " = ", kind, "(", paste(parts, collapse = ", "), ")", sep = "")
}

convert_method <- function(e, name_of_factors)
{
  key = names(e)[1]
  r = e[[1]]
  settings = c(mcmc_termination = "termination", smc_termination = "termination", adaptive_target = "target",
               adaptive_resampling = "resampling", enk_shifter = "shifter")
  if (key %in% names(settings))
  {
    b = convert_builtin_call(r)
    return(list(where = "algorithm", line = paste("  ", settings[[key]], " = ", convert_call_text(b$name, b$arguments, r), sep = "")))
  }
  if (key %in% c("smc_sequence", "enk_sequence"))
  {
    b = convert_builtin_call(r)
    arguments = sapply(b$arguments, convert_argument, record = r)
    # the third argument: the factors affected, by number
    if (length(b$arguments) == 3 && !grepl("=", b$arguments[3]))
    {
      numbers = eval(parse(text = b$arguments[3]))
      arguments[3] = paste("factors =", name_of_factors(numbers, r))
    }
    return(list(where = "algorithm", line = paste("  sequence = ", b$name, "(", paste(arguments, collapse = ", "), ")", sep = "")))
  }
  if (key == "mcmc_weights")
  {
    if (!isTRUE(r$custom) && grepl("^c\\([-+0-9., eE]*\\)$", gsub("\\s+", "", r$block_function)))
      return(list(where = "algorithm", line = paste("  weights =", convert_argument(r$block_function, r))))
    return(list(where = "algorithm", line = paste("  weights =", convert_function(r))))
  }
  if (key == "ssm")
  {
    b = convert_builtin_call(r)
    return(list(where = "model", line = paste("  time =", convert_call_text(b$name, b$arguments, r))))
  }
  convert_stop(r, paste(key, "blocks are not in the ilike language yet."))
}
