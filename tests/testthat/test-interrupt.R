# Long runs can be interrupted (dev/plan.md 3.11): the C++ code asks R for
# user interrupts between chunks, between SMC steps and within long chains.
# A time limit is raised by the same check, so it stands in for a user's
# interrupt. Rcpp turns it into an R interrupt condition, which testthat would
# take as an interrupt of the tests, so each case runs in its own R process.

run_interrupted <- function(call_text)
{
  script = tempfile(fileext = ".R")
  writeLines(c(
    "suppressMessages(library(ilike))",
    "ext = function(f) system.file('extdata', f, package = 'ilike')",
    "algorithm = tempfile(fileext = '.ilike')",
    "writeLines(c('/***m_proposal,ilike::norm_rw(theta,1)***/', '', '/***mcmc_termination,ilike::iterations(50000000)***/'), algorithm)",
    "mcmc = suppressMessages(ilike::compile(c(ext('test_gaussian_posterior.ilike'), algorithm)))",
    "is = suppressMessages(ilike::compile(ext('test_gaussian_posterior.ilike')))",
    "smc = suppressMessages(ilike::compile(c(ext('test_gaussian_posterior.ilike'), ext('test_smc_adaptive_cess.ilike'))))",
    "f = function() { on.exit(setTimeLimit(), add = TRUE); setTimeLimit(elapsed = 1, transient = TRUE); suppressMessages(", call_text, ") }",
    "t = system.time(r <- tryCatch({ f(); 'finished' }, interrupt = function(c) 'interrupted', error = function(e) conditionMessage(e)))[['elapsed']]",
    "cat('RESULT', r, t, '\\n')",
    "cat('AFTER', suppressMessages(ilike::IS(is, number_of_importance_points = 1000, seed = 2)), '\\n')"), script)
  output = system2(file.path(R.home("bin"), "Rscript"), script, stdout = TRUE, stderr = TRUE)
  output
}

test_that("long MCMC, parallel IS and SMC runs stop at an interrupt, and the next run works", {
  skip_on_cran() # slow (over 3 s); runs on CI and locally
  # under AddressSanitizer the child R processes would start without the
  # sanitizer runtime (Rscript drops DYLD_INSERT_LIBRARIES)
  skip_if(nzchar(Sys.getenv("ASAN_OPTIONS")), "not under AddressSanitizer")

  calls = c("ilike::MCMC(mcmc, initial_values = list(list(theta = 1)), seed = 1)",
            "ilike::IS(is, number_of_importance_points = 1e8, seed = 1, parallel = TRUE)",
            "ilike::SMC_with_MCMC(smc, number_of_particles = 5e6, seed = 1)")
  for (call_text in calls)
  {
    output = run_interrupted(call_text)
    result = strsplit(grep("^RESULT", output, value = TRUE), " ")[[1]]
    expect_true(result[2] %in% c("interrupted", "reached"), label = call_text)
    expect_lt(as.numeric(result[length(result)]), 10)
    after = as.numeric(strsplit(grep("^AFTER", output, value = TRUE), " ")[[1]][2])
    expect_equal(after, dnorm(1, 1, sqrt(2), log = TRUE), tolerance = 0.05)
  }
})
