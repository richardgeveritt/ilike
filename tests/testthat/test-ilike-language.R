# The parser of the ilike language (dev/ilike-language.md, I3a): the syntax
# tree of a file, and errors with file:line:column.

language_tree = function(text) ilike:::ilike_language_syntax_tree(text, "test.ilike", FALSE)

test_that("the ilike language is recognised by its version line, and old files are not", {
  expect_true(ilike:::is_ilike_language_text("// a comment\n/* another */\nilike 1.0\nmodel { }"))
  old = paste(readLines(system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")), collapse = "\n")
  expect_false(ilike:::is_ilike_language_text(old))
  expect_false(ilike:::is_ilike_language_text("model { }"))
})

test_that("sections, statements and expressions are read into the syntax tree", {
  text = paste(c(
    "ilike 1.0",
    "import \"ssm.ilike\" as inner",
    "constants {",
    "  matrix Sigma0",
    "  real prior_sd = 10   // a default",
    "  vector[2] m0 = [0, 0]",
    "}",
    "parameters { real<lower=0> sigma; any tree }",
    "model {",
    "  theta ~ norm(0, prior_sd)",
    "  prior p = density(log = cpp::log_prior,",
    "                    simulate = cpp::simulate_prior)",
    "  factor lik = cpp::log_likelihood",
    "}",
    "cpp {",
    "  // a brace in a comment } and in strings \"}\" '}'",
    "  double f(const Parameters &parameters) { return 1.0; }",
    "}",
    "algorithm mcmc_rw {",
    "  move m1 = metropolis(mvnorm_rw(beta, cov = Sigma0), adapt = scaled_am(averaged))",
    "  sweep = deterministic(m1); termination = iterations(2.38^2 / 2)",
    "}"), collapse = "\n")
  tree = language_tree(text)
  expect_equal(tree$version, "1.0")
  expect_equal(sapply(tree$sections, function(s) s$section),
               c("import", "constants", "parameters", "model", "cpp", "algorithm"))
  expect_equal(tree$sections[[1]]$alias, "inner")

  constants = tree$sections[[2]]$statements
  expect_equal(sapply(constants, function(s) s$keyword), c("matrix", "real", "vector"))
  expect_equal(constants[[2]]$value$value, 10)
  expect_equal(constants[[3]]$type_arguments[[1]]$value, 2)
  expect_equal(constants[[3]]$value$expression, "vector")

  parameters = tree$sections[[3]]$statements
  expect_equal(parameters[[1]]$constraints$lower$value, 0)
  expect_equal(parameters[[2]]$keyword, "any")

  model = tree$sections[[4]]$statements
  expect_equal(sapply(model, function(s) s$operator), c("~", "=", "="))
  expect_equal(sapply(model, function(s) s$keyword), c("", "prior", "factor"))
  expect_equal(model[[1]]$value$`function`, "norm")
  expect_equal(names(model[[2]]$value$arguments), c("log", "simulate"))
  expect_equal(model[[2]]$value$arguments$log$namespace, "cpp")
  expect_equal(model[[3]]$line, 13)

  expect_match(tree$sections[[5]]$code, "double f\\(const Parameters &parameters\\) \\{ return 1.0; \\}")
  expect_match(tree$sections[[5]]$code, "'\\}'")

  algorithm = tree$sections[[6]]
  expect_equal(algorithm$name, "mcmc_rw")
  expect_equal(length(algorithm$statements), 3)
  expect_equal(names(algorithm$statements[[1]]$value$arguments), c("", "adapt"))
  expect_equal(algorithm$statements[[3]]$value$arguments[[1]]$operator, "/")
})

test_that("loops and deterministic statements are parsed (I4)", {
  tree = language_tree("ilike 1.0\nmodel {\n  for (i in 1:n) {\n    mu[i] <- a + b * x[i]\n    y[i] ~ norm(mu[i], 1)\n  }\n  z<-1\n}")
  statements = tree$sections[[1]]$statements
  expect_equal(statements[[1]]$keyword, "for")
  expect_equal(statements[[1]]$index, "i")
  expect_equal(statements[[1]]$range$operator, ":")
  expect_equal(length(statements[[1]]$body), 2)
  expect_equal(statements[[1]]$body[[1]]$operator, "<-")
  expect_equal(statements[[1]]$body[[2]]$operator, "~")
  # "<-" only when the two characters are adjacent
  expect_equal(statements[[2]]$operator, "<-")
  expect_error(language_tree("ilike 1.0\nmodel {\n  w < -1\n}"), "test.ilike:3:7: expected a constraint")
})

test_that("terms of a named factor are parsed", {
  tree = language_tree("ilike 1.0\nmodel {\n  for (e in 1:n) factor ising += b * s[e]\n  factor one += 2\n}")
  statements = tree$sections[[1]]$statements
  term = statements[[1]]$body[[1]]
  expect_equal(term$keyword, "factor")
  expect_equal(term$target$name, "ising")
  expect_equal(term$operator, "+=")
  expect_equal(term$value$operator, "*")
  expect_equal(statements[[2]]$operator, "+=")
})

test_that("syntax errors give the file, line and column", {
  expect_error(language_tree("ilike 1.0\nmodel { theta ~ norm(0, 1 }"), "test.ilike:2:27: expected \",\" or \"\\)\"")
  expect_error(language_tree("ilike 1.0\nmodle { }"), "test.ilike:2:1: unknown section \"modle\"")
  expect_error(language_tree("ilike 2.0\n"), "reads version 1 of the ilike language")
  expect_error(language_tree("model { }"), "begins with the version line")
  expect_error(language_tree("ilike 1.0\nmodel {\n  for (i on 1:n) y[i] ~ norm(0, 1)\n}"), "test.ilike:3:10: expected \"in\" after the loop's index")
  expect_error(language_tree("ilike 1.0\nmodel {\n  for (i in 1:n) {\n y[i] ~ norm(0, 1)\n"), "test.ilike:3:18: the loop's block has no closing brace")
  expect_error(language_tree("ilike 1.0\ncpp { int f() { return 1; }\n"), "the cpp section has no closing brace")
  expect_error(language_tree("ilike 1.0\nmodel { theta norm(0, 1) }"), "test.ilike:2:21: expected \"=\" or \"~\"|expected")
})
