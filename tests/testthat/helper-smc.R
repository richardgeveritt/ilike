test_smc_z <- function(algorithm_file, results_name, seed = 1,
                       model_file = "test_gaussian_posterior.ilike")
{
  # SMC with MCMC moves on model_file, with the SMC settings in algorithm_file
  filename_model = system.file("extdata", model_file, package = "ilike")
  filename_algorithm = system.file("extdata", algorithm_file, package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_algorithm))
  ilike::SMC_with_MCMC(recipe, results_name = results_name, number_of_particles = 2000,
                       smc_iterations_to_store = 100, seed = seed)
}

test_smc_annealing_z <- function(results_name, seed = 1)
{
  # fixed annealing schedule c(0,0.25,0.5,0.75,1)
  test_smc_z("test_smc_annealing.ilike", results_name, seed)
}

# log evidence of prior N(1,1) times the N(1 | theta,1) likelihood raised to power p
annealed_log_evidence <- function(p)
{
  -0.5*p*log(2*pi) - 0.5*log(1+p)
}

# Weighted mean and variance of one scalar parameter at each target, using
# each particle's last stored MCMC iteration (the state the weight belongs to).
smc_weighted_moments_by_target <- function(output)
{
  t(sapply(sort(unique(output$Target)), function(target) {
    o = output[output$Target == target, ]
    o = o[o$Iteration == max(o$Iteration), ]
    w = exp(o$LogWeight - max(o$LogWeight))
    w = w/sum(w)
    m = sum(w*o$Value)
    c(target = target, mean = m, var = sum(w*(o$Value - m)^2))
  }))
}
