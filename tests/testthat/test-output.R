# Wide and long output formats (dev/plan.md S7): load_*_output(wide = TRUE)
# gives one row per point; ilike_pivot_longer turns it back into the long
# format (one row per scalar) that the load functions give by default.

expect_same_output <- function(a, b)
{
  a = as.data.frame(a); b = as.data.frame(b)
  expect_identical(names(a), names(b))
  expect_equal(a, b, tolerance = 1e-12, ignore_attr = TRUE)
}

test_that("wide output and ilike_pivot_longer round-trip for MCMC, SMC and PF", {
  on.exit(unlink(c("output_mcmc", "output_smc", "output_pf"), recursive = TRUE), add = TRUE)
  ext = function(f) system.file("extdata", f, package = "ilike")
  pl = list(c(1, 1), diag(2))
  r = ilike::compile(c(ext("test_gaussian_posterior_vector.ilike"), ext("test_imh_proposal_mvnorm_2d.ilike")), model_parameter_list = pl)
  ilike::MCMC(r, model_parameter_list = pl, results_name = "output_mcmc", seed = 1,
              initial_values = list(list(theta = c(1, 1)), list(theta = c(0, 0))))
  r = ilike::compile(c(ext("test_gaussian_posterior.ilike"), ext("test_smc_annealing.ilike")))
  ilike::SMC_with_MCMC(r, results_name = "output_smc", number_of_particles = 50, smc_iterations_to_store = 100, seed = 1)
  ilike::PF(ilike::compile(ext("test_ssm_one_line.ilike")), results_name = "output_pf", number_of_particles = 40, seed = 1)

  long = ilike::load_mcmc_output("output_mcmc"); wide = ilike::load_mcmc_output("output_mcmc", wide = TRUE)
  expect_equal(nrow(wide), nrow(long)/2)
  expect_true(all(c("theta_1", "theta_2", "Chain", "Iteration") %in% names(wide)))
  expect_same_output(wide, ilike::ilike_pivot_wider(long))
  expect_same_output(ilike::ilike_pivot_longer(wide), long)

  for (d in c("output_smc", "output_pf"))
  {
    long = ilike::load_smc_output(d); wide = ilike::load_smc_output(d, wide = TRUE)
    expect_same_output(wide, ilike::ilike_pivot_wider(long))
    expect_same_output(ilike::ilike_pivot_longer(wide), long)
  }
  expect_error(ilike::load_smc_output("output_smc", wide = TRUE, ggmcmc = TRUE, ggsmc = FALSE, ilike.output = FALSE, as.mcmc = TRUE), "incompatible")
})
