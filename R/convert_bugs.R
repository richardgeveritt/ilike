# Converting BUGS, JAGS and NIMBLE models to the ilike language (dev/plan.md
# 3.4, C231). BUGS code is R-like: after the model { } wrapper and
# truncation T(a, b) are rewritten, R's parser reads it (NIMBLE code is R
# already), and the statements are walked as R expressions. Each dialect's
# parameterisations are mapped to ilike's (precisions to standard
# deviations, BUGS's argument orders, NIMBLE's named alternatives).

#' Convert a BUGS, JAGS or NIMBLE model to the ilike language
#'
#' Translates a model written in the BUGS language (as used by WinBUGS,
#' OpenBUGS and JAGS) or in NIMBLE's dialect of it into a model in the ilike
#' language, with a data section holding the data given (which also tells
#' which variables are data: the others are parameters). Translated: stochastic
#' relations (~), deterministic relations (<-, with logit(p) <- … and
#' log(x) <- … on the left), for loops, arithmetic and the functions exp,
#' log, sqrt, pow, abs, logit, ilogit (and NIMBLE's expit), sin, cos, tanh,
#' min, max, equals, step, and inprod and sum over vectors whose lengths are known from
#' the data; truncation T(a, b) (and BUGS's I(a, b)) for dnorm. Distributions
#' (with BUGS's parameterisations, or NIMBLE's when dialect = "nimble"):
#' dnorm, dlnorm, dt, dgamma, dbeta, dunif, dexp, dchisqr (dchisq), dweib,
#' ddexp, dlogis, dpois, dbern, dbin, dnegbin, dcat, dmulti, ddirch (ddirich),
#' dmnorm, dmt and dwish (whose precision matrices must be data: they are
#' inverted into new data, name_inverse). Missing values in the data,
#' censoring (dinterval), other truncations, precision matrices that are
#' parameters, and other functions stop with a message.
#'
#' An algorithm section is added with a NUTS move on all the parameters when
#' they are all continuous and unconstrained vectors or scalars.
#'
#' @param model The model: a file name, the code as a character string, or an
#'   R expression (for example from nimbleCode() or quote()).
#' @param data A list with the data (and constants), by name.
#' @param dialect (optional) "jags" (the default), "bugs" or "nimble": the
#'   parameterisations of the distributions (JAGS and BUGS share them).
#' @param output (optional) The file to write.
#' @param overwrite (optional) If TRUE, an existing output file is replaced.
#' @return Invisibly, the converted text.
#' @export
convert_bugs_model <- function(model, data = list(), dialect = c("jags", "bugs", "nimble"), output = NULL, overwrite = FALSE)
{
  dialect = match.arg(dialect)
  if (is.language(model))
    code = model
  else
  {
    if ((length(model) == 1) && file.exists(model))
      text = paste(readLines(model, warn = FALSE), collapse = "\n")
    else
      text = paste(model, collapse = "\n")
    code = bugs_parse(text)
  }
  if (is.call(code) && identical(code[[1]], as.name("quote"))) code = code[[2]]
  lines = bugs_translate(code, data, dialect)
  if (!is.null(output))
  {
    if (file.exists(output) && !overwrite)
      stop("convert_bugs_model: ", output, " exists (use overwrite = TRUE to replace it).")
    writeLines(lines, output)
  }
  invisible(paste(lines, collapse = "\n"))
}

bugs_parse <- function(text)
{
  text = gsub("#[^\n]*", "", text)
  if (grepl("(^|[^A-Za-z0-9_.])data[[:space:]]*\\{", text))
    stop("convert_bugs_model: data { } blocks (data transformations in JAGS) cannot be translated: transform the data in R.")
  text = sub("^[[:space:]]*(model|var)[[:space:]]*\\{", "{", text)
  # truncation: dist(...) T(a, b) -> dist(...) %T% c(a, b)
  text = gsub("\\)[[:space:]]*[TI][[:space:]]*\\(", ") %T% c(", text)
  parsed = tryCatch(parse(text = text, keep.source = FALSE), error = function(e)
    stop("convert_bugs_model: the model could not be read: ", conditionMessage(e)))
  if (length(parsed) != 1) stop("convert_bugs_model: the model is one block model { ... }.")
  parsed[[1]]
}

bugs_missing <- function(x) identical(x, quote(expr = ))
# which of a call's index arguments are empty (x[, j]) or ranges (x[1:K]),
# by position (an empty argument cannot be passed to a function)
bugs_empty_at <- function(indices) vapply(seq_along(indices), function(j) identical(indices[[j]], quote(expr = )), logical(1))
bugs_range_at <- function(indices) vapply(seq_along(indices), function(j) {
  if (identical(indices[[j]], quote(expr = ))) return(FALSE)
  i = indices[[j]]; is.call(i) && identical(i[[1]], as.name(":")) }, logical(1))

bugs_translate <- function(code, data, dialect)
{
  tr = new.env()
  tr$data = data
  tr$extra_data = list()
  tr$dialect = dialect
  tr$statements = character(0)
  tr$depth = 1
  tr$parameters = character(0)
  tr$discrete = FALSE
  tr$constrained = FALSE
  tr$nodes = character(0)
  for (name in names(data))
    if (any(is.na(data[[name]])))
      stop("convert_bugs_model: the data ", name, " has missing values (partly observed variables cannot be translated).")
  bugs_statement(tr, code)
  all_data = c(data, tr$extra_data)
  items = vapply(names(all_data), function(name) paste0(name, " = ", bugs_data_value(all_data[[name]], name)), character(1))
  out = c("ilike 1.0", paste0("// converted from ", toupper(dialect), " by convert_bugs_model()"))
  if (length(items) > 0)
    out = c(out, paste0("data { ", paste(items, collapse = "; "), " }"))
  out = c(out, "model {", tr$statements, "}")
  parameters = unique(tr$parameters)
  if ((length(parameters) > 0) && !tr$discrete && !tr$constrained)
    out = c(out, "algorithm mcmc {",
            paste0("  move m = nuts([", paste(parameters, collapse = ", "), "], step = 0.1, adapt = dual_averaging(until = iterations(1000)))"),
            "  termination = iterations(3000)", "}")
  else if (length(parameters) > 0)
    out = c(out, "// (no algorithm section: some parameters are discrete, which nuts cannot move)")
  out
}

bugs_data_value <- function(x, name)
{
  if (is.logical(x)) x = as.numeric(x)
  if (!is.numeric(x)) stop("convert_bugs_model: the data ", name, " is not numeric.")
  number = function(v) vapply(v, function(u) { t = format(u, digits = 15); if (as.numeric(t) == u) t else sprintf("%.17g", u) }, character(1))
  if (is.matrix(x))
    return(paste0("[", paste(apply(x, 1, function(r) paste0("[", paste(number(r), collapse = ", "), "]")), collapse = ", "), "]"))
  if (!is.null(dim(x)) && length(dim(x)) > 2)
    stop("convert_bugs_model: the data ", name, " has more than two dimensions.")
  if (length(x) == 1) return(number(x))
  paste0("[", paste(number(x), collapse = ", "), "]")
}

bugs_line <- function(tr, text) tr$statements = c(tr$statements, paste0(strrep("  ", tr$depth), text))

bugs_statement <- function(tr, e)
{
  if (!is.call(e)) stop("convert_bugs_model: unexpected statement ", deparse(e), ".")
  head = as.character(e[[1]])
  if (head == "{")
  {
    for (k in seq_along(e)[-1]) bugs_statement(tr, e[[k]])
    return(invisible())
  }
  if (head == "for")
  {
    range = e[[3]]
    if (!is.call(range) || as.character(range[[1]]) != ":") stop("convert_bugs_model: a loop's range is a:b.")
    bugs_line(tr, paste0("for (", as.character(e[[2]]), " in ", bugs_expr(tr, range[[2]]), ":", bugs_expr(tr, range[[3]]), ") {"))
    tr$depth = tr$depth + 1
    bugs_statement(tr, e[[4]])
    tr$depth = tr$depth - 1
    bugs_line(tr, "}")
    return(invisible())
  }
  if (head == "~")
    return(bugs_stochastic(tr, e[[2]], e[[3]]))
  if (head %in% c("<-", "="))
    return(bugs_deterministic(tr, e[[2]], e[[3]]))
  stop("convert_bugs_model: the statement ", paste(deparse(e), collapse = " "), " cannot be translated.")
}

bugs_name <- function(e)
{
  if (is.name(e)) return(as.character(e))
  if (is.call(e) && identical(e[[1]], as.name("["))) return(bugs_name(e[[2]]))
  NULL
}

# a whole variable written x, x[] or x[, ] or x[1:K] (for vector arguments)
bugs_whole <- function(e)
{
  if (is.name(e)) return(as.character(e))
  if (is.call(e) && identical(e[[1]], as.name("[")) && is.name(e[[2]]))
  {
    indices = as.list(e)[-(1:2)]
    if (all(bugs_empty_at(indices) | bugs_range_at(indices)))
      return(as.character(e[[2]]))
  }
  NULL
}

bugs_deterministic <- function(tr, lhs, rhs)
{
  value = bugs_expr(tr, rhs)
  if (is.call(lhs) && !identical(lhs[[1]], as.name("[")))
  {
    link = as.character(lhs[[1]])
    if (link == "logit") value = paste0("inv_logit(", value, ")")
    else if (link == "log") value = paste0("exp(", value, ")")
    else stop("convert_bugs_model: the link function ", link, " on the left of <- cannot be translated (logit and log can).")
    lhs = lhs[[2]]
  }
  name = bugs_name(lhs)
  if (!is.null(tr$data[[name]])) stop("convert_bugs_model: ", name, " is data and also defined by <- (transform the data in R).")
  tr$nodes = c(tr$nodes, name)
  bugs_line(tr, paste0(bugs_expr(tr, lhs), " <- ", value))
}

bugs_functions <- c(exp = "exp", log = "log", sqrt = "sqrt", pow = "pow", abs = "abs", logit = "logit", ilogit = "inv_logit",
                    expit = "inv_logit", sin = "sin", cos = "cos", tanh = "tanh", min = "min", max = "max",
                    log1p = "log1p", expm1 = "expm1", lgamma = "lgamma", loggam = "lgamma", equals = "equals", step = "step")

# the number of elements of a data variable's dimension k
bugs_extent <- function(tr, name, k)
{
  x = tr$data[[name]]
  if (is.null(x)) return(NULL)
  if (is.matrix(x)) return(dim(x)[k])
  if (k == 1) return(length(x))
  NULL
}

bugs_number <- function(tr, e)
{
  if (is.numeric(e)) return(e)
  if (is.name(e) && !is.null(tr$data[[as.character(e)]]) && length(tr$data[[as.character(e)]]) == 1) return(as.numeric(tr$data[[as.character(e)]]))
  NULL
}

# the elements of a vector argument of sum/inprod: x[] / x[1:K] / X[i, ] /
# X[i, 1:K], as a function of k
bugs_vector_terms <- function(tr, e, known = NULL)
{
  if (is.name(e)) e = call("[", e, quote(expr = ))
  if (!is.call(e) || !identical(e[[1]], as.name("["))) stop("convert_bugs_model: sum and inprod take vectors written x[] or x[a:b].")
  name = as.character(e[[2]])
  indices = as.list(e)[-(1:2)]
  open = which(bugs_empty_at(indices) | bugs_range_at(indices))
  if (length(open) != 1) stop("convert_bugs_model: sum and inprod take one open index (x[], X[i, ], x[1:K]).")
  if (bugs_empty_at(indices)[open])
  {
    from = 1
    to = bugs_extent(tr, name, if (length(indices) == 2) open else 1)
  }
  else
  {
    i = indices[[open]]
    from = bugs_number(tr, i[[2]]); to = bugs_number(tr, i[[3]])
  }
  if (is.null(to) && !is.null(known) && !is.null(from)) to = from + known - 1
  if (is.null(from) || is.null(to)) stop("convert_bugs_model: the length of ", name, " in sum or inprod must be known from the data.")
  ks = seq(from, to)
  vapply(ks, function(k) { idx = indices; idx[[open]] = k
    paste0(name, "[", paste(vapply(seq_along(idx), function(j) bugs_expr(tr, idx[[j]]), character(1)), collapse = ", "), "]") }, character(1))
}

bugs_expr <- function(tr, e)
{
  paren = function(x) paste0("(", x, ")")
  if (is.numeric(e))
  {
    t = format(e, digits = 15)
    return(if (as.numeric(t) == e) t else sprintf("%.17g", e))
  }
  if (is.name(e)) return(as.character(e))
  if (!is.call(e)) stop("convert_bugs_model: unexpected ", deparse(e), ".")
  head = as.character(e[[1]])
  if (head == "(") return(paren(bugs_expr(tr, e[[2]])))
  if (head == "[")
  {
    indices = as.list(e)[-(1:2)]
    if (any(bugs_empty_at(indices)))
      stop("convert_bugs_model: ", paste(deparse(e), collapse = ""), ": a whole row or column is translated in sum, inprod and vector arguments only.")
    return(paste0(bugs_expr(tr, e[[2]]), "[", paste(vapply(seq_along(indices), function(j) bugs_expr(tr, indices[[j]]), character(1)), collapse = ", "), "]"))
  }
  if (head %in% c("+", "-", "*", "/", "^"))
  {
    if (length(e) == 2) return(paste0(head, paren(bugs_expr(tr, e[[2]]))))
    return(paste0(paren(bugs_expr(tr, e[[2]])), " ", head, " ", paren(bugs_expr(tr, e[[3]]))))
  }
  if (head == "sum")
    return(paren(paste(bugs_vector_terms(tr, e[[2]]), collapse = " + ")))
  if (head == "inprod")
  {
    # (the length from whichever vector is data)
    a = tryCatch(bugs_vector_terms(tr, e[[2]]), error = function(err) NULL)
    b = bugs_vector_terms(tr, e[[3]], if (is.null(a)) NULL else length(a))
    if (is.null(a)) a = bugs_vector_terms(tr, e[[2]], length(b))
    if (length(a) != length(b)) stop("convert_bugs_model: inprod of vectors of different lengths.")
    return(paren(paste(paste0(a, " * ", b), collapse = " + ")))
  }
  f = bugs_functions[head]
  if (is.na(f)) stop("convert_bugs_model: the function ", head, " cannot be translated.")
  paste0(f, "(", paste(vapply(as.list(e)[-1], function(a) bugs_expr(tr, a), character(1)), collapse = ", "), ")")
}

# a precision matrix that is data, as the name of its inverse (new data)
bugs_inverse <- function(tr, e, what)
{
  name = bugs_whole(e)
  if (is.null(name) || is.null(tr$data[[name]]) || !is.matrix(tr$data[[name]]))
    stop("convert_bugs_model: the ", what, " matrix of ", paste(deparse(e), collapse = ""),
         " must be data (it is inverted when converting: ilike's statements take covariance matrices).")
  inverse = paste0(name, "_inverse")
  tr$extra_data[[inverse]] = solve(tr$data[[name]])
  inverse
}

bugs_stochastic <- function(tr, lhs, rhs)
{
  truncation = NULL
  if (is.call(rhs) && identical(rhs[[1]], as.name("%T%")))
  {
    bounds = rhs[[3]]
    truncation = list(lower = if (length(bounds) >= 2 && !bugs_missing(bounds[[2]])) bugs_expr(tr, bounds[[2]]) else "-Inf",
                      upper = if (length(bounds) >= 3 && !bugs_missing(bounds[[3]])) bugs_expr(tr, bounds[[3]]) else "Inf")
    rhs = rhs[[2]]
  }
  if (!is.call(rhs)) stop("convert_bugs_model: a distribution is needed after ~.")
  name = bugs_name(lhs)
  observed = !is.null(tr$data[[name]])
  if (!observed) tr$parameters = c(tr$parameters, name)
  distribution = as.character(rhs[[1]])
  arguments = as.list(rhs)[-1]
  keys = names(arguments)
  if (is.null(keys)) keys = rep("", length(arguments))
  named = function(key) if (key %in% keys) arguments[[which(keys == key)]] else NULL
  positional = arguments[keys == ""]
  pos = function(k) if (length(positional) >= k) positional[[k]] else NULL
  x = function(e) bugs_expr(tr, e)
  whole_arg = function(e, what) {
    w = bugs_whole(e)
    if (is.null(w)) stop("convert_bugs_model: the ", what, " of ", distribution, " is a variable written whole (x or x[]).")
    w
  }
  nimble = (tr$dialect == "nimble")
  inv_sqrt = function(e) if (is.numeric(e)) bugs_expr(tr, 1/sqrt(e)) else paste0("1 / sqrt(", x(e), ")")
  inv = function(e) if (is.numeric(e)) bugs_expr(tr, 1/e) else paste0("1 / (", x(e), ")")
  # a scale-type argument: given as sd/scale, var, or precision/rate
  sd_of = function(position, sd_key = "sd", var_key = "var", tau_key = "tau") {
    if (!is.null(named(sd_key))) return(x(named(sd_key)))
    if (!is.null(named(var_key))) return(paste0("sqrt(", x(named(var_key)), ")"))
    if (!is.null(named(tau_key))) return(inv_sqrt(named(tau_key)))
    inv_sqrt(pos(position))
  }
  whole = FALSE
  call = switch(distribution,
    dnorm = paste0("norm(", x(if (is.null(named("mean"))) pos(1) else named("mean")), ", ", sd_of(2), ")"),
    dlnorm = paste0("lnorm(", x(if (is.null(named("meanlog"))) pos(1) else named("meanlog")), ", ", sd_of(2, "sdlog", "varlog", "taulog"), ")"),
    dt = paste0("t(", x(if (is.null(named("df"))) pos(3) else named("df")), ", ", x(if (is.null(named("mu"))) pos(1) else named("mu")), ", ",
                if (!is.null(named("sigma"))) x(named("sigma")) else inv_sqrt(if (is.null(named("tau"))) pos(2) else named("tau")), ")"),
    dgamma = paste0("gamma(", x(if (is.null(named("shape"))) pos(1) else named("shape")), ", ",
                    if (!is.null(named("scale"))) inv(named("scale")) else x(if (is.null(named("rate"))) pos(2) else named("rate")), ")"),
    dbeta = paste0("beta(", x(pos(1)), ", ", x(pos(2)), ")"),
    dunif = paste0("unif(", x(pos(1)), ", ", x(pos(2)), ")"),
    dexp = paste0("exp(", if (!is.null(named("scale"))) inv(named("scale")) else x(if (is.null(named("rate"))) pos(1) else named("rate")), ")"),
    dchisqr = , dchisq = paste0("chisq(", x(pos(1)), ")"),
    dweib = if (nimble) paste0("weibull(", x(pos(1)), ", ",
                                if (!is.null(named("rate"))) inv(named("rate")) else if (!is.null(named("lambda"))) paste0("pow(", x(named("lambda")), ", -1 / (", x(pos(1)), "))")
                                else x(if (is.null(named("scale"))) pos(2) else named("scale")), ")")
            else paste0("weibull(", x(pos(1)), ", pow(", x(pos(2)), ", -1 / (", x(pos(1)), ")))"),
    ddexp = paste0("laplace(", x(pos(1)), ", ", if (nimble) (if (!is.null(named("rate"))) inv(named("rate")) else x(if (is.null(named("scale"))) pos(2) else named("scale"))) else inv(pos(2)), ")"),
    dlogis = paste0("logis(", x(pos(1)), ", ", if (nimble) (if (!is.null(named("rate"))) inv(named("rate")) else x(if (is.null(named("scale"))) pos(2) else named("scale"))) else inv(pos(2)), ")"),
    dpois = paste0("pois(", x(if (is.null(named("lambda"))) pos(1) else named("lambda")), ")"),
    dbern = paste0("bern(", x(if (is.null(named("prob"))) pos(1) else named("prob")), ")"),
    dbin = , dbinom = paste0("binom(", x(if (is.null(named("size"))) pos(2) else named("size")), ", ", x(if (is.null(named("prob"))) pos(1) else named("prob")), ")"),
    dnegbin = , dnbinom = paste0("negbinom(", x(if (is.null(named("size"))) pos(2) else named("size")), ", ", x(if (is.null(named("prob"))) pos(1) else named("prob")), ")"),
    dcat = paste0("categorical(", whole_arg(if (is.null(named("prob"))) pos(1) else named("prob"), "probabilities"), ")"),
    dmulti = , dmultinom = { whole = TRUE; paste0("multinom(", x(if (is.null(named("size"))) pos(2) else named("size")), ", ", whole_arg(if (is.null(named("prob"))) pos(1) else named("prob"), "probabilities"), ")") },
    ddirch = , ddirich = { whole = TRUE; paste0("dirichlet(", whole_arg(if (is.null(named("alpha"))) pos(1) else named("alpha"), "concentrations"), ")") },
    dmnorm = { whole = TRUE
      covariance = if (!is.null(named("cov"))) whole_arg(named("cov"), "covariance")
                   else bugs_inverse(tr, if (!is.null(named("prec"))) named("prec") else pos(2), "precision")
      paste0("mvnorm(", whole_arg(if (is.null(named("mean"))) pos(1) else named("mean"), "mean"), ", ", covariance, ")") },
    dmt = { whole = TRUE
      paste0("mvt(", x(if (is.null(named("df"))) pos(3) else named("df")), ", ", whole_arg(if (is.null(named("mu"))) pos(1) else named("mu"), "location"), ", ",
             if (!is.null(named("scale"))) whole_arg(named("scale"), "scale") else bugs_inverse(tr, if (!is.null(named("prec"))) named("prec") else pos(2), "precision"), ")") },
    dwish = { whole = TRUE
      paste0("wishart(", x(if (is.null(named("df"))) pos(2) else named("df")), ", ",
             if (!is.null(named("S"))) whole_arg(named("S"), "scale") else bugs_inverse(tr, if (!is.null(named("R"))) named("R") else pos(1), "scale"), ")") },
    dinterval = stop("convert_bugs_model: censoring (dinterval) cannot be translated."),
    stop("convert_bugs_model: the distribution ", distribution, " cannot be translated."))
  if (!is.null(truncation))
  {
    if (distribution != "dnorm") stop("convert_bugs_model: truncation T(,) is translated for dnorm only (", name, ").")
    call = paste0("tnorm(", sub("^norm\\(", "", sub("\\)$", "", call)), ", ", truncation$lower, ", ", truncation$upper, ")")
  }
  if (!observed)
  {
    if (distribution %in% c("dpois", "dbern", "dbin", "dbinom", "dnegbin", "dnbinom", "dcat", "dmulti", "dmultinom")) tr$discrete = TRUE
    # (Dirichlet and Wishart parameters: nuts moves them through its automatic transforms, C232)
  }
  target = if (whole) { w = bugs_whole(lhs); if (is.null(w)) stop("convert_bugs_model: the target of ", distribution, " is a variable written whole."); w } else bugs_expr(tr, lhs)
  bugs_line(tr, paste0(target, " ~ ", call))
}
