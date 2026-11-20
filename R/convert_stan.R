# Converting Stan programs to the ilike language (dev/plan.md 3.4, C230).
# A tokenizer and a recursive-descent parser for the subset of Stan described
# in convert_stan_model's documentation, and a translation of its data,
# transformed data, parameters, transformed parameters and model blocks into
# a data section and model statements. Vectorised statements become loops
# (their lengths from the declarations), matrix-vector products become sums
# (their lengths from the data), and the constraints of parameters become
# the supports of their priors.

#' Convert a Stan program to the ilike language
#'
#' Translates a Stan program into a model in the ilike language, with a data
#' section holding the data given, so that it can be run by any of ilike's
#' algorithms. What is translated: the blocks data, transformed data,
#' parameters, transformed parameters and model (generated quantities are left
#' out, with a warning; functions blocks are refused); declarations of int,
#' real, vector, row_vector, matrix, simplex, cov_matrix and corr_matrix
#' variables, arrays of them (array[N] real and the older real y[N]), with
#' lower and upper bounds; for loops; assignments (as deterministic nodes);
#' sampling statements, vectorised ones included (as loops over the elements),
#' with truncation T[a, b] for normal; target += with a log density function
#' (y ~ dist(...) written as target += dist_lpdf(y | ...)) or an expression;
#' arithmetic, element-wise functions (exp, log, sqrt, square, inv_logit,
#' logit, inv, fabs, abs, log1p, expm1, lgamma, sin, cos, tanh, fmin, fmax,
#' pow), and matrix-vector products X * beta whose inner length is known from
#' the data. Distributions: normal, std_normal, lognormal, gamma, inv_gamma, beta,
#' exponential, uniform, cauchy, student_t, double_exponential, logistic,
#' weibull, chi_square, von_mises, poisson, poisson_log, bernoulli,
#' bernoulli_logit, binomial, binomial_logit, neg_binomial, neg_binomial_2,
#' categorical, multinomial, dirichlet, multi_normal, multi_student_t,
#' wishart, inv_wishart, lkj_corr. A constrained parameter's prior becomes the
#' distribution on its support (normal with bounds: tnorm; normal(0, s),
#' cauchy(0, s) and student_t(nu, 0, s) with lower = 0: halfnorm, halfcauchy,
#' halft); a parameter with both bounds and no prior gets the uniform
#' distribution, a simplex the uniform Dirichlet and a correlation matrix
#' lkj(1). Anything else (if, while, user-defined functions, improper priors,
#' other constraints) stops with a message naming the statement.
#'
#' The ilike model has proper priors and normalised densities, so its
#' evidence is that of the model; Stan's target drops constants, which does
#' not change the posterior. An algorithm section is added with a NUTS move on
#' all the parameters (with dual averaging of the step over the first 1000
#' iterations), when they can all be moved by it.
#'
#' @param model The Stan program: a file name, or the program as a character string.
#' @param data A list with the data, named as in the data block.
#' @param output (optional) The file to write.
#' @param overwrite (optional) If TRUE, an existing output file is replaced.
#' @return Invisibly, the converted text.
#' @export
convert_stan_model <- function(model, data = list(), output = NULL, overwrite = FALSE)
{
  if ((length(model) == 1) && file.exists(model))
    text = paste(readLines(model, warn = FALSE), collapse = "\n")
  else
    text = paste(model, collapse = "\n")
  tokens = stan_tokens(text)
  program = stan_parse_program(tokens)
  lines = stan_translate(program, data)
  if (!is.null(output))
  {
    if (file.exists(output) && !overwrite)
      stop("convert_stan_model: ", output, " exists (use overwrite = TRUE to replace it).")
    writeLines(lines, output)
  }
  invisible(paste(lines, collapse = "\n"))
}

# --- tokens ---------------------------------------------------------------------

stan_tokens <- function(text)
{
  # comments: /* ... */, // ... and # ...
  text = gsub("/\\*.*?\\*/", " ", text, perl = TRUE)
  text = gsub("//[^\n]*", "", text)
  text = gsub("#[^\n]*", "", text)
  pattern = paste0("(?:[0-9]+\\.?[0-9]*(?:[eE][-+]?[0-9]+)?|\\.[0-9]+(?:[eE][-+]?[0-9]+)?)",
                   "|[A-Za-z_][A-Za-z0-9_]*",
                   "|<=|>=|==|!=|&&|\\|\\||\\+=|-=|\\*=|/=|\\.\\*|\\./|<-",
                   "|[-+*/^%<>=!?:;,()\\[\\]{}~|']")
  matches = regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]]
  rest = gsub(pattern, "", text, perl = TRUE)
  if (grepl("[^[:space:]]", rest))
    stop("convert_stan_model: unexpected characters in the program: ", substr(gsub("[[:space:]]", "", rest), 1, 20))
  matches
}

# a cursor over the tokens
stan_cursor <- function(tokens)
{
  state = new.env()
  state$tokens = tokens
  state$i = 1
  state
}
stan_peek <- function(s, k = 0) if (s$i + k <= length(s$tokens)) s$tokens[s$i + k] else ""
stan_next <- function(s) { t = stan_peek(s); s$i = s$i + 1; t }
stan_expect <- function(s, t)
{
  got = stan_next(s)
  if (got != t)
    stop("convert_stan_model: expected \"", t, "\" but found \"", got, "\" near: ",
         paste(s$tokens[max(1, s$i - 6):min(length(s$tokens), s$i + 3)], collapse = " "))
}
stan_accept <- function(s, t) { if (stan_peek(s) == t) { s$i = s$i + 1; TRUE } else FALSE }
stan_is_name <- function(t) grepl("^[A-Za-z_][A-Za-z0-9_]*$", t)
stan_is_number <- function(t) grepl("^([0-9]|\\.[0-9])", t)

# --- parser ---------------------------------------------------------------------

stan_parse_program <- function(tokens)
{
  s = stan_cursor(tokens)
  blocks = list()
  while (s$i <= length(tokens))
  {
    name = stan_next(s)
    if ((name == "transformed") || (name == "generated"))
      name = paste(name, stan_next(s))
    if (name == "functions")
      stop("convert_stan_model: functions blocks (user-defined functions) cannot be translated.")
    if (!(name %in% c("data", "transformed data", "parameters", "transformed parameters", "model", "generated quantities")))
      stop("convert_stan_model: unknown block \"", name, "\".")
    stan_expect(s, "{")
    statements = list()
    while (stan_peek(s) != "}")
      statements[[length(statements) + 1]] = stan_parse_statement(s)
    stan_expect(s, "}")
    blocks[[name]] = statements
  }
  blocks
}

stan_types <- c("int", "real", "vector", "row_vector", "matrix", "simplex", "unit_vector", "ordered", "positive_ordered",
                "cov_matrix", "corr_matrix", "cholesky_factor_cov", "cholesky_factor_corr", "array", "complex")

stan_parse_statement <- function(s)
{
  t = stan_peek(s)
  if (t == "{")
  {
    stan_next(s)
    body = list()
    while (stan_peek(s) != "}")
      body[[length(body) + 1]] = stan_parse_statement(s)
    stan_expect(s, "}")
    return(list(kind = "block", body = body))
  }
  if (t == "for")
  {
    stan_next(s); stan_expect(s, "(")
    variable = stan_next(s); stan_expect(s, "in")
    from = stan_parse_expression(s); stan_expect(s, ":"); to = stan_parse_expression(s)
    stan_expect(s, ")")
    return(list(kind = "for", variable = variable, from = from, to = to, body = stan_parse_statement(s)))
  }
  if (t %in% c("if", "while", "break", "continue", "return", "reject", "fatal_error"))
    stop("convert_stan_model: \"", t, "\" statements cannot be translated.")
  if (t == "print")
  {
    depth = 0
    repeat { x = stan_next(s); if (x == "(") depth = depth + 1; if (x == ")") depth = depth - 1; if ((x == ";") && (depth == 0)) break }
    return(list(kind = "skip"))
  }
  if (t %in% stan_types)
    return(stan_parse_declaration(s))
  if (t == "target")
  {
    stan_next(s); stan_expect(s, "+=")
    value = stan_parse_expression(s); stan_expect(s, ";")
    return(list(kind = "target", value = value))
  }
  lhs = stan_parse_expression(s)
  op = stan_next(s)
  if (op == "~")
  {
    distribution = stan_next(s)
    stan_expect(s, "(")
    arguments = stan_parse_arguments(s, ")")
    truncation = NULL
    if (stan_peek(s) == "T" && stan_peek(s, 1) == "[")
    {
      stan_next(s); stan_next(s)
      lower = if (stan_peek(s) == ",") NULL else stan_parse_expression(s)
      stan_expect(s, ",")
      upper = if (stan_peek(s) == "]") NULL else stan_parse_expression(s)
      stan_expect(s, "]")
      truncation = list(lower = lower, upper = upper)
    }
    stan_expect(s, ";")
    return(list(kind = "sample", target = lhs, distribution = distribution, arguments = arguments, truncation = truncation))
  }
  if (op %in% c("=", "<-"))
  {
    value = stan_parse_expression(s); stan_expect(s, ";")
    return(list(kind = "assign", target = lhs, value = value))
  }
  stop("convert_stan_model: the statement with \"", op, "\" cannot be translated (only =, ~ and target +=).")
}

stan_parse_arguments <- function(s, close)
{
  arguments = list()
  if (stan_accept(s, close))
    return(arguments)
  repeat
  {
    arguments[[length(arguments) + 1]] = stan_parse_expression(s)
    if (stan_accept(s, "|")) { arguments[[length(arguments) + 1]] = list(kind = "bar"); next }
    if (stan_accept(s, close)) break
    stan_expect(s, ",")
  }
  arguments
}

# type<constraints>[sizes] name[dims] (= value);  or  array[dims] type<...>[sizes] name;
stan_parse_declaration <- function(s)
{
  array_dims = list()
  type = stan_next(s)
  if (type == "array")
  {
    stan_expect(s, "[")
    array_dims = stan_parse_arguments(s, "]")
    type = stan_next(s)
  }
  if (type %in% c("unit_vector", "ordered", "positive_ordered", "cholesky_factor_cov", "cholesky_factor_corr", "complex"))
    stop("convert_stan_model: ", type, " variables cannot be translated.")
  lower = NULL; upper = NULL
  if (stan_accept(s, "<"))
  {
    repeat
    {
      key = stan_next(s); stan_expect(s, "=")
      value = stan_parse_expression(s, no_greater = TRUE)
      if (key == "lower") lower = value
      else if (key == "upper") upper = value
      else if (!(key %in% c("offset", "multiplier"))) stop("convert_stan_model: unknown constraint ", key, ".")
      if (stan_accept(s, ">")) break
      stan_expect(s, ",")
    }
  }
  sizes = list()
  if (stan_accept(s, "["))
    sizes = stan_parse_arguments(s, "]")
  name = stan_next(s)
  if (stan_accept(s, "["))
    array_dims = c(array_dims, stan_parse_arguments(s, "]"))
  value = NULL
  if (stan_accept(s, "=") || stan_accept(s, "<-"))
    value = stan_parse_expression(s)
  stan_expect(s, ";")
  list(kind = "declare", name = name, type = type, lower = lower, upper = upper, sizes = sizes, array_dims = array_dims, value = value)
}

# expressions, by precedence
stan_parse_expression <- function(s, no_greater = FALSE)
{
  e = stan_parse_binary(s, 1, no_greater)
  if (stan_accept(s, "?"))
  {
    a = stan_parse_expression(s); stan_expect(s, ":"); b = stan_parse_expression(s)
    stop("convert_stan_model: conditional expressions (a ? b : c) cannot be translated.")
  }
  e
}

stan_levels <- list(c("||"), c("&&"), c("==", "!="), c("<", "<=", ">", ">="), c("+", "-"), c("*", "/", "%", ".*", "./"))

stan_parse_binary <- function(s, level, no_greater)
{
  if (level > length(stan_levels))
    return(stan_parse_unary(s, no_greater))
  left = stan_parse_binary(s, level + 1, no_greater)
  repeat
  {
    op = stan_peek(s)
    if (!(op %in% stan_levels[[level]]) || (no_greater && (op == ">")))
      break
    stan_next(s)
    right = stan_parse_binary(s, level + 1, no_greater)
    left = list(kind = "binary", op = op, left = left, right = right)
  }
  left
}

stan_parse_unary <- function(s, no_greater)
{
  if (stan_accept(s, "-")) return(list(kind = "unary", op = "-", operand = stan_parse_unary(s, no_greater)))
  if (stan_accept(s, "+")) return(stan_parse_unary(s, no_greater))
  if (stan_accept(s, "!")) stop("convert_stan_model: logical negation cannot be translated.")
  base = stan_parse_postfix(s)
  if (stan_accept(s, "^"))
    return(list(kind = "binary", op = "^", left = base, right = stan_parse_unary(s, no_greater)))
  base
}

stan_parse_postfix <- function(s)
{
  e = stan_parse_primary(s)
  repeat
  {
    if (stan_peek(s) == "[")
    {
      stan_next(s)
      indices = list()
      repeat
      {
        if (stan_peek(s) %in% c(",", "]")) indices[[length(indices) + 1]] = list(kind = "all")
        else
        {
          a = stan_parse_expression(s)
          if (stan_accept(s, ":")) a = list(kind = "range", from = a, to = stan_parse_expression(s))
          indices[[length(indices) + 1]] = a
        }
        if (stan_accept(s, "]")) break
        stan_expect(s, ",")
      }
      e = list(kind = "index", target = e, indices = indices)
    }
    else if (stan_accept(s, "'"))
      e = list(kind = "transpose", operand = e)
    else
      break
  }
  e
}

stan_parse_primary <- function(s)
{
  t = stan_next(s)
  if (stan_is_number(t))
    return(list(kind = "number", value = as.numeric(t), text = t))
  if (t == "(")
  {
    e = stan_parse_expression(s); stan_expect(s, ")")
    return(e)
  }
  if (t == "[")
    return(list(kind = "vector", elements = stan_parse_arguments(s, "]")))
  if (t == "{")
    return(list(kind = "vector", elements = stan_parse_arguments(s, "}")))
  if (stan_is_name(t))
  {
    if (stan_accept(s, "("))
      return(list(kind = "call", name = t, arguments = stan_parse_arguments(s, ")")))
    return(list(kind = "name", name = t))
  }
  stop("convert_stan_model: unexpected \"", t, "\".")
}

# --- translation ----------------------------------------------------------------

stan_translate <- function(program, data)
{
  if (!is.null(program[["generated quantities"]]))
    warning("convert_stan_model: the generated quantities block is left out.")
  tr = new.env()
  tr$data = data
  tr$declarations = list()   # name -> declaration, with role
  tr$statements = character(0)
  tr$depth = 1
  tr$indices = character(0)  # loop indices in scope
  tr$counter = 0
  tr$priors = character(0)   # parameters with a prior statement
  tr$factor_count = 0
  for (st in program[["data"]])
  {
    if (st$kind != "declare") stop("convert_stan_model: the data block holds declarations only.")
    st$role = "data"
    tr$declarations[[st$name]] = st
    if (is.null(data[[st$name]]))
      stop("convert_stan_model: the data ", st$name, " is not in the data given.")
  }
  stan_emit(tr, "// transformed data", program[["transformed data"]], "transformed data")
  for (st in program[["parameters"]])
  {
    if (st$kind != "declare") stop("convert_stan_model: the parameters block holds declarations only.")
    if (!is.null(st$value)) stop("convert_stan_model: a parameter cannot be assigned (", st$name, ").")
    st$role = "parameter"
    tr$declarations[[st$name]] = st
  }
  stan_emit(tr, "// transformed parameters", program[["transformed parameters"]], "transformed parameters")
  stan_emit(tr, "// model", program[["model"]], "model")
  # parameters without a prior statement
  defaults = character(0)
  for (name in names(tr$declarations))
  {
    d = tr$declarations[[name]]
    if ((d$role != "parameter") || (name %in% tr$priors))
      next
    defaults = c(defaults, stan_default_prior(tr, d))
  }
  # the data section
  data_names = names(tr$declarations)[vapply(tr$declarations, function(d) d$role == "data", logical(1))]
  data_items = vapply(data_names, function(name) paste0(name, " = ", stan_data_value(data[[name]], name)), character(1))
  parameters = names(tr$declarations)[vapply(tr$declarations, function(d) d$role == "parameter", logical(1))]
  out = c("ilike 1.0", "// converted from Stan by convert_stan_model()")
  if (length(data_items) > 0)
    out = c(out, paste0("data { ", paste(data_items, collapse = "; "), " }"))
  out = c(out, "model {", if (length(defaults) > 0) c("  // priors from the parameters' constraints (none in the program)", defaults),
          tr$statements, "}")
  # (simplex, cov_matrix and corr_matrix through nuts's automatic transforms, C232)
  movable = all(vapply(parameters, function(name) tr$declarations[[name]]$type %in%
                         c("real", "vector", "row_vector", "matrix", "simplex", "cov_matrix", "corr_matrix"), logical(1)))
  if (movable && (length(parameters) > 0))
    out = c(out, "algorithm mcmc {",
            paste0("  move m = nuts([", paste(parameters, collapse = ", "), "], step = 0.1, adapt = dual_averaging(until = iterations(1000)))"),
            "  termination = iterations(3000)", "}")
  else if (length(parameters) > 0)
    out = c(out, "// (no algorithm section: some parameters are of types nuts cannot move)")
  out
}

stan_data_value <- function(x, name)
{
  if (is.logical(x)) x = as.numeric(x)
  if (!is.numeric(x)) stop("convert_stan_model: the data ", name, " is not numeric.")
  if (any(is.na(x))) stop("convert_stan_model: the data ", name, " has missing values.")
  # the shortest text that reads back as the same number
  number = function(v) vapply(v, function(u) { t = format(u, digits = 15); if (as.numeric(t) == u) t else sprintf("%.17g", u) }, character(1))
  if (is.matrix(x))
    return(paste0("[", paste(apply(x, 1, function(r) paste0("[", paste(number(r), collapse = ", "), "]")), collapse = ", "), "]"))
  if (!is.null(dim(x)) && length(dim(x)) > 2)
    stop("convert_stan_model: the data ", name, " has more than two dimensions.")
  if (length(x) == 1) return(number(x))
  paste0("[", paste(number(x), collapse = ", "), "]")
}

# the numeric value of a size expression, from the data
stan_value <- function(tr, e)
{
  if (e$kind == "number") return(e$value)
  if (e$kind == "name" && !is.null(tr$data[[e$name]]) && length(tr$data[[e$name]]) == 1) return(as.numeric(tr$data[[e$name]]))
  if (e$kind == "binary" && e$op %in% c("+", "-", "*", "/"))
  {
    a = stan_value(tr, e$left); b = stan_value(tr, e$right)
    if (is.null(a) || is.null(b)) return(NULL)
    return(switch(e$op, "+" = a + b, "-" = a - b, "*" = a*b, "/" = a/b))
  }
  NULL
}

stan_line <- function(tr, text) tr$statements = c(tr$statements, paste0(strrep("  ", tr$depth), text))

stan_new_index <- function(tr)
{
  repeat
  {
    tr$counter = tr$counter + 1
    name = paste0("i", tr$counter)
    if (is.null(tr$declarations[[name]]) && !(name %in% tr$indices)) return(name)
  }
}

stan_emit <- function(tr, comment, statements, block)
{
  if (length(statements) == 0) return(invisible())
  stan_line(tr, comment)
  for (st in statements) stan_statement(tr, st, block)
}

stan_statement <- function(tr, st, block)
{
  switch(st$kind,
    skip = invisible(),
    block = for (b in st$body) stan_statement(tr, b, block),
    declare = {
      if (block %in% c("data", "parameters"))
        stop("convert_stan_model: unexpected declaration of ", st$name, ".")
      st$role = if (block == "transformed data") "data_node" else "node"
      tr$declarations[[st$name]] = st
      if (!is.null(st$value))
        stan_assign(tr, list(kind = "name", name = st$name), st$value)
    },
    assign = stan_assign(tr, st$target, st$value),
    "for" = {
      if (!is.null(tr$declarations[[st$variable]])) stop("convert_stan_model: the loop index ", st$variable, " is also a variable.")
      stan_line(tr, paste0("for (", st$variable, " in ", stan_expr(tr, st$from), ":", stan_expr(tr, st$to), ") {"))
      tr$indices = c(tr$indices, st$variable)
      tr$depth = tr$depth + 1
      stan_statement(tr, st$body, block)
      tr$depth = tr$depth - 1
      tr$indices = setdiff(tr$indices, st$variable)
      stan_line(tr, "}")
    },
    sample = {
      if (block != "model") stop("convert_stan_model: sampling statements belong in the model block.")
      stan_sample(tr, st$target, st$distribution, st$arguments, st$truncation)
    },
    target = {
      if (block != "model") stop("convert_stan_model: target += belongs in the model block.")
      v = st$value
      if ((v$kind == "call") && grepl("_(lpdf|lpmf|log)$", v$name))
      {
        bar = which(vapply(v$arguments, function(a) a$kind == "bar", logical(1)))
        if (length(bar) != 1) stop("convert_stan_model: target += ", v$name, "(y | ...) needs its | .")
        distribution = sub("_(lpdf|lpmf|log)$", "", v$name)
        stan_sample(tr, v$arguments[[1]], distribution, v$arguments[-(1:bar)], NULL)
      }
      else
      {
        if (stan_length(tr, v) != "") stop("convert_stan_model: target += with a vector expression cannot be translated (it needs a sum).")
        tr$factor_count = tr$factor_count + 1
        stan_line(tr, paste0("factor target", tr$factor_count, " += ", stan_expr(tr, v)))
      }
    },
    stop("convert_stan_model: a statement of kind ", st$kind, " cannot be translated."))
}

# the declaration of the variable an expression refers to (name or indexed name)
stan_base_name <- function(e)
{
  if (e$kind == "name") return(e$name)
  if (e$kind == "index") return(stan_base_name(e$target))
  NULL
}

# dimensions of a declared variable: array dims then vector/matrix sizes
stan_dims <- function(d)
{
  sizes = d$sizes
  if (d$type %in% c("simplex", "cov_matrix", "corr_matrix")) sizes = if (d$type == "simplex") d$sizes else c(d$sizes, d$sizes)
  c(d$array_dims, sizes)
}

# the length (as ilike text) of a one-dimensional expression; "" for a
# scalar; stops for a matrix-valued expression used element-wise
stan_length <- function(tr, e)
{
  switch(e$kind,
    number = "",
    name = {
      if (e$name %in% tr$indices) return("")
      d = tr$declarations[[e$name]]
      if (is.null(d)) stop("convert_stan_model: unknown variable ", e$name, ".")
      dims = stan_dims(d)
      if (length(dims) == 0) return("")
      if (length(dims) == 1) return(stan_expr(tr, dims[[1]]))
      return("matrix")
    },
    index = {
      d = tr$declarations[[stan_base_name(e)]]
      dims = stan_dims(d)
      left = length(dims) - length(e$indices)
      ranges = vapply(e$indices, function(i) i$kind %in% c("range", "all"), logical(1))
      if (any(ranges)) stop("convert_stan_model: slices (x[a:b], x[, j]) cannot be translated.")
      if (left <= 0) return("")
      if (left == 1) return(stan_expr(tr, dims[[length(dims)]]))
      "matrix"
    },
    unary = stan_length(tr, e$operand),
    binary = {
      if ((e$op == "*") && (stan_length(tr, e$left) == "matrix"))
      {
        # a matrix times a vector: a vector with the matrix's rows
        d = tr$declarations[[stan_base_name(e$left)]]
        return(stan_expr(tr, stan_dims(d)[[1]]))
      }
      a = stan_length(tr, e$left); b = stan_length(tr, e$right)
      if (a == "matrix" || b == "matrix") stop("convert_stan_model: matrix arithmetic other than a matrix times a vector cannot be translated.")
      if (a != "") a else b
    },
    call = {
      if (is.na(stan_functions[e$name]) && !(e$name %in% c("square", "inv")))
        stop("convert_stan_model: the function ", e$name, " cannot be translated.")
      ls = vapply(e$arguments, function(a) stan_length(tr, a), character(1))
      ls = ls[ls != ""]
      if (length(ls) == 0) "" else ls[1]
    },
    transpose = stan_length(tr, e$operand),
    stop("convert_stan_model: an expression of kind ", e$kind, " cannot be translated."))
}

stan_functions <- c(exp = "exp", log = "log", sqrt = "sqrt", inv_logit = "inv_logit", logit = "logit", fabs = "abs", abs = "abs",
                    log1p = "log1p", expm1 = "expm1", lgamma = "lgamma", sin = "sin", cos = "cos", tanh = "tanh",
                    fmin = "min", fmax = "max", pow = "pow")

# an expression as ilike text; `element` (a loop index) takes element
# `element` of one-dimensional parts (for vectorised statements)
stan_expr <- function(tr, e, element = NULL)
{
  paren = function(x) paste0("(", x, ")")
  switch(e$kind,
    number = e$text,
    name = {
      if (!is.null(element) && stan_length(tr, e) != "") paste0(e$name, "[", element, "]") else e$name
    },
    index = {
      indices = vapply(e$indices, function(i) stan_expr(tr, i), character(1))
      if (!is.null(element) && stan_length(tr, e) != "") indices = c(indices, element)
      paste0(stan_expr(tr, e$target), "[", paste(indices, collapse = ", "), "]")
    },
    unary = if (e$operand$kind == "number") paste0("-", e$operand$text) else paste0("-", paren(stan_expr(tr, e$operand, element))),
    binary = {
      if ((e$op == "*") && (stan_length(tr, e$left) == "matrix"))
      {
        if (is.null(element)) stop("convert_stan_model: a matrix-vector product is translated element by element only (in a vectorised statement or an assignment to a vector).")
        X = e$left
        d = tr$declarations[[stan_base_name(X)]]
        K = stan_value(tr, stan_dims(d)[[2]])
        if (is.null(K)) stop("convert_stan_model: the number of columns of ", stan_base_name(X), " must be known from the data.")
        rowname = stan_expr(tr, X)
        terms = vapply(seq_len(K), function(k)
          paste0(rowname, "[", element, ", ", k, "] * ", stan_expr(tr, list(kind = "index", target = e$right, indices = list(list(kind = "number", value = k, text = as.character(k)))))),
          character(1))
        return(paren(paste(terms, collapse = " + ")))
      }
      op = switch(e$op, ".*" = "*", "./" = "/", e$op)
      if (op %in% c("||", "&&", "==", "!=", "<", "<=", ">", ">=", "%"))
        stop("convert_stan_model: the operator ", e$op, " cannot be translated.")
      paste0(paren(stan_expr(tr, e$left, element)), " ", op, " ", paren(stan_expr(tr, e$right, element)))
    },
    call = {
      args = vapply(e$arguments, function(a) stan_expr(tr, a, element), character(1))
      if (e$name == "square") return(paste0("pow(", args[1], ", 2)"))
      if (e$name == "inv") return(paste0("1 / (", args[1], ")"))
      if (e$name %in% c("sum", "mean", "dot_product", "rep_vector", "rep_array", "segment", "head", "tail", "cumulative_sum"))
        stop("convert_stan_model: the function ", e$name, " cannot be translated.")
      f = stan_functions[e$name]
      if (is.na(f)) stop("convert_stan_model: the function ", e$name, " cannot be translated.")
      paste0(f, "(", paste(args, collapse = ", "), ")")
    },
    vector = paste0("[", paste(vapply(e$elements, function(a) stan_expr(tr, a), character(1)), collapse = ", "), "]"),
    stop("convert_stan_model: an expression of kind ", e$kind, " cannot be translated."))
}

stan_assign <- function(tr, target, value)
{
  name = stan_base_name(target)
  d = tr$declarations[[name]]
  if (is.null(d) || !(d$role %in% c("node", "data_node")))
    stop("convert_stan_model: only local and transformed variables can be assigned (", name, ").")
  n = stan_length(tr, target)
  if (n == "matrix") stop("convert_stan_model: assigning a whole matrix (", name, ") cannot be translated.")
  if (n == "")
  {
    stan_line(tr, paste0(stan_expr(tr, target), " <- ", stan_expr(tr, value)))
    return(invisible())
  }
  i = stan_new_index(tr)
  stan_line(tr, paste0("for (", i, " in 1:", n, ") ", stan_expr(tr, target, i), " <- ", stan_expr(tr, value, i)))
}

# Stan's distributions: the ilike call (as text) from the arguments' texts
stan_distribution <- function(name, a, n_trials = NULL)
{
  f = function(g, ...) paste0(g, "(", paste(c(...), collapse = ", "), ")")
  switch(name,
    normal = f("norm", a[1], a[2]), std_normal = "norm(0, 1)", lognormal = f("lnorm", a[1], a[2]), gamma = f("gamma", a[1], a[2]),
    inv_gamma = f("invgamma", a[1], a[2]), beta = f("beta", a[1], a[2]), exponential = f("exp", a[1]),
    uniform = f("unif", a[1], a[2]), cauchy = f("cauchy", a[1], a[2]), student_t = f("t", a[1], a[2], a[3]),
    double_exponential = f("laplace", a[1], a[2]), logistic = f("logis", a[1], a[2]), weibull = f("weibull", a[1], a[2]),
    chi_square = f("chisq", a[1]), von_mises = f("vonmises", a[1], a[2]),
    poisson = f("pois", a[1]), poisson_log = f("pois", paste0("exp(", a[1], ")")),
    bernoulli = f("bern", a[1]), bernoulli_logit = f("bern", paste0("inv_logit(", a[1], ")")),
    binomial = f("binom", a[1], a[2]), binomial_logit = f("binom", a[1], paste0("inv_logit(", a[2], ")")),
    neg_binomial = f("negbinom", a[1], paste0("(", a[2], ") / (1 + (", a[2], "))")),
    neg_binomial_2 = f("negbinom", a[2], paste0("(", a[2], ") / ((", a[1], ") + (", a[2], "))")),
    categorical = f("categorical", a[1]), multinomial = f("multinom", n_trials, a[1]), dirichlet = f("dirichlet", a[1]),
    multi_normal = f("mvnorm", a[1], a[2]), multi_student_t = f("mvt", a[1], a[2], a[3]),
    wishart = f("wishart", a[1], a[2]), inv_wishart = f("invwishart", a[1], a[2]),
    lkj_corr = f("lkj", a[1]),
    stop("convert_stan_model: the distribution ", name, " cannot be translated."))
}

stan_vector_distributions <- c("categorical", "multinomial", "dirichlet", "multi_normal", "multi_student_t", "wishart", "inv_wishart", "lkj_corr")

stan_sample <- function(tr, target, distribution, arguments, truncation)
{
  name = stan_base_name(target)
  d = tr$declarations[[name]]
  if (is.null(d)) stop("convert_stan_model: unknown variable ", name, ".")
  if (d$role %in% c("node", "data_node"))
    stop("convert_stan_model: a sampling statement on the transformed variable ", name, " (a Jacobian would be needed) cannot be translated.")
  if (d$role == "parameter") tr$priors = c(tr$priors, name)
  whole = distribution %in% stan_vector_distributions && distribution != "categorical"
  n = if (whole) "" else stan_length(tr, target)
  if (n == "matrix") stop("convert_stan_model: a vectorised statement on the matrix ", name, " cannot be translated.")
  element = if (n == "") NULL else stan_new_index(tr)
  args = vapply(arguments, function(a) {
    if (whole || distribution == "categorical")
    {
      if (!(a$kind %in% c("name", "number", "vector"))) stop("convert_stan_model: the vector arguments of ", distribution, " must be variables or numbers (", name, ").")
      return(stan_expr(tr, a))
    }
    stan_expr(tr, a, element)
  }, character(1))
  n_trials = NULL
  if (distribution == "multinomial")
  {
    if (is.null(tr$data[[name]])) stop("convert_stan_model: multinomial needs its counts (", name, ") as data.")
    n_trials = as.character(sum(tr$data[[name]]))
  }
  if (distribution == "lkj_corr")
  {
    K = stan_value(tr, stan_dims(d)[[1]])
    args = c(args, if (is.null(K)) character(0) else as.character(K))
  }
  call = stan_distribution(distribution, args, n_trials)
  if (distribution == "lkj_corr" && length(args) == 2) call = paste0("lkj(", args[1], ", ", args[2], ")")
  # truncation, and the parameter's constraints
  lower = if (!is.null(truncation$lower)) stan_expr(tr, truncation$lower) else NULL
  upper = if (!is.null(truncation$upper)) stan_expr(tr, truncation$upper) else NULL
  if (!is.null(truncation) && distribution != "normal")
    stop("convert_stan_model: truncation T[,] is translated for normal only (", name, ").")
  if ((d$role == "parameter") && is.null(truncation))
  {
    lower = if (!is.null(d$lower)) stan_expr(tr, d$lower) else NULL
    upper = if (!is.null(d$upper)) stan_expr(tr, d$upper) else NULL
  }
  if (!is.null(lower) || !is.null(upper))
    call = stan_constrained(name, distribution, args, lower, upper, call)
  targettext = if (is.null(element)) stan_expr(tr, target) else stan_expr(tr, target, element)
  if (is.null(element))
    stan_line(tr, paste0(targettext, " ~ ", call))
  else
    stan_line(tr, paste0("for (", element, " in 1:", n, ") ", targettext, " ~ ", call))
}

# a distribution on a variable with bounds: the distribution on that support
stan_constrained <- function(name, distribution, a, lower, upper, call)
{
  zero = function(x) !is.null(x) && suppressWarnings(!is.na(as.numeric(x))) && as.numeric(x) == 0
  one = function(x) !is.null(x) && suppressWarnings(!is.na(as.numeric(x))) && as.numeric(x) == 1
  positive = c("lognormal", "gamma", "inv_gamma", "exponential", "weibull", "chi_square")
  if ((distribution %in% positive) && zero(lower) && is.null(upper)) return(call)
  if ((distribution == "beta") && zero(lower) && one(upper)) return(call)
  if ((distribution == "uniform") && !is.null(lower) && !is.null(upper) && (lower == a[1]) && (upper == a[2])) return(call)
  if (distribution %in% c("dirichlet", "wishart", "inv_wishart", "lkj_corr")) return(call)
  if (distribution == "std_normal")
  {
    distribution = "normal"
    a = c("0", "1")
  }
  if (distribution == "normal")
  {
    if (zero(lower) && is.null(upper) && zero(a[1])) return(paste0("halfnorm(", a[2], ")"))
    return(paste0("tnorm(", a[1], ", ", a[2], ", ", if (is.null(lower)) "-Inf" else lower, ", ", if (is.null(upper)) "Inf" else upper, ")"))
  }
  if ((distribution == "cauchy") && zero(lower) && is.null(upper) && zero(a[1])) return(paste0("halfcauchy(", a[2], ")"))
  if ((distribution == "student_t") && zero(lower) && is.null(upper) && zero(a[2])) return(paste0("halft(", a[1], ", ", a[3], ")"))
  stop("convert_stan_model: the prior of ", name, " is ", distribution, " on a constrained range (lower = ", if (is.null(lower)) "none" else lower,
       ", upper = ", if (is.null(upper)) "none" else upper, "), which has no counterpart in ilike (it would need a truncated ", distribution, ").")
}

# a parameter with no sampling statement: the prior its constraints imply, if proper
stan_default_prior <- function(tr, d)
{
  name = d$name
  dims = stan_dims(d)
  if (d$type == "simplex")
  {
    K = stan_value(tr, d$sizes[[1]])
    if (is.null(K)) stop("convert_stan_model: the size of the simplex ", name, " must be known from the data.")
    return(paste0("  ", name, " ~ dirichlet([", paste(rep("1", K), collapse = ", "), "])"))
  }
  if (d$type == "corr_matrix")
  {
    K = stan_value(tr, d$sizes[[1]])
    return(paste0("  ", name, " ~ lkj(1", if (is.null(K)) "" else paste0(", ", K), ")"))
  }
  if (!is.null(d$lower) && !is.null(d$upper) && (length(dims) <= 1))
  {
    call = paste0("unif(", stan_expr(tr, d$lower), ", ", stan_expr(tr, d$upper), ")")
    if (length(dims) == 0) return(paste0("  ", name, " ~ ", call))
    return(paste0("  for (i0 in 1:", stan_expr(tr, dims[[1]]), ") ", name, "[i0] ~ ", call))
  }
  stop("convert_stan_model: the parameter ", name, " has no prior (Stan's improper flat prior has no counterpart in ilike): give it one in the model block.")
}
