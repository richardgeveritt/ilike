# Long runs can be interrupted (dev/plan.md 3.11): the C++ code asks R for
# user interrupts between chunks, between SMC steps and within long chains.
# Rcpp turns an interrupt into an R interrupt condition, which testthat would
# take as an interrupt of the tests, so each case runs in its own R process.
#
# On Unix the child process sends itself a real SIGINT (what Ctrl-C does)
# after one second. Windows has no such signal, so there an elapsed time limit
# stands in for it; R raises a time limit only once the run has checked
# several times after the limit passed, so on Windows only the chain is
# tested (it checks every few thousand iterations; a parallel IS with one
# long wave checks too rarely for the time limit, though not for SIGINT).

run_interrupted <- function(call_text)
{
  if (.Platform$OS.type == "windows")
    start = "setTimeLimit(elapsed = 1, transient = TRUE)"
  else
    # the whole of 'sleep 1; kill' goes to the background, so the signal
    # arrives while the run is going (system() ignores SIGINT while it waits)
    start = "system2('sh', c('-c', shQuote(paste('sleep 1; kill -INT', Sys.getpid()))), wait = FALSE)"
  script = tempfile(fileext = ".R")
  writeLines(c(
    "suppressMessages(library(ilike))",
    "ext = function(f) system.file('extdata', f, package = 'ilike')",
    "algorithm = tempfile(fileext = '.ilike')",
    "writeLines(c('/***m_proposal,ilike::norm_rw(theta,1)***/', '', '/***mcmc_termination,ilike::iterations(50000000)***/'), algorithm)",
    "mcmc = suppressMessages(ilike::compile(c(ext('test_gaussian_posterior.ilike'), algorithm)))",
    "is = suppressMessages(ilike::compile(ext('test_gaussian_posterior.ilike')))",
    "smc_algorithm = tempfile(fileext = '.ilike')",
    "writeLines(c('/***m_proposal,ilike::norm_rw(theta,1)***/', '', '/***mcmc_termination,ilike::iterations(2000)***/', '', '/***smc_sequence,ilike::annealing(power,c(0,1))***/', '', '/***adaptive_target,ilike::cess(0.95)***/'), smc_algorithm)",
    "smc = suppressMessages(ilike::compile(c(ext('test_gaussian_posterior.ilike'), smc_algorithm)))",
    "nested = suppressMessages(ilike::compile(ext('test_algorithmic_pf.ilike')))",
    paste0("f = function() { on.exit(setTimeLimit(), add = TRUE); ", start, "; suppressMessages(", call_text, ") }"),
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

  # long runs (each takes well over a second) that use little memory:
  # storage for every point is reserved at the start, so a huge IS or an
  # unthinned chain would spend its time allocating, and page on a small CI
  # machine
  calls = c("ilike::MCMC(mcmc, initial_values = list(list(theta = 1)), seed = 1, thin = 1000)",
            "ilike::IS(nested, number_of_importance_points = 50000, seed = 1, parallel = TRUE)",
            "ilike::SMC_with_MCMC(smc, number_of_particles = 20000, seed = 1)")
  if (.Platform$OS.type == "windows")
    calls = calls[1]
  for (call_text in calls)
  {
    output = run_interrupted(call_text)
    result = strsplit(grep("^RESULT", output, value = TRUE), " ")[[1]]
    expect_true(result[2] %in% c("interrupted", "reached"), label = call_text)
    expect_lt(as.numeric(result[length(result)]), 30)
    after = as.numeric(strsplit(grep("^AFTER", output, value = TRUE), " ")[[1]][2])
    expect_equal(after, dnorm(1, 1, sqrt(2), log = TRUE), tolerance = 0.05)
  }
})
