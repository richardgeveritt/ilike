# Test custom Metropolis proposal in MCMC

run_mcmc_m_proposal <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_m_proposal",seed=2,initial_values = list(list(theta=1)))
}

# Test custom Metropolis proposal in MCMC where only the first factor (the prior) is used

run_mcmc_m_proposal_prior <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal_prior.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_m_proposal_prior",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis proposal in MCMC
# ilike::norm_rw

run_mcmc_m_proposal_norm_rw <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal_norm_rw.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_m_proposal_norm_rw",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis proposal in MCMC
# ilike::unif_rw

run_mcmc_m_proposal_unif_rw <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal_unif_rw.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_m_proposal_unif_rw",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis proposal in MCMC
# ilike::mirror

run_mcmc_m_proposal_mirror <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal_mirror.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_m_proposal_mirror",seed=2,initial_values = list(list(theta=1)))
}

# Test custom independent Metropolis-Hastings proposal in MCMC

run_mcmc_imh_proposal <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_imh_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_imh_proposal",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in independent Metropolis-Hastings proposal in MCMC
# ilike::norm

run_mcmc_imh_proposal_norm <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_imh_proposal_norm.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_imh_proposal_norm",seed=2,initial_values = list(list(theta=1)))
}

# Test custom Metropolis-Hastings proposal in MCMC

run_mcmc_mh_proposal_custom <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_custom",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis-Hastings proposal in MCMC
# ilike::norm_rw

run_mcmc_mh_proposal_norm_rw <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal_norm_rw.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_norm_rw",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis-Hastings proposal in MCMC
# ilike::unif_rw

run_mcmc_mh_proposal_unif_rw <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal_unif_rw.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_unif_rw",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis-Hastings proposal in MCMC
# ilike::langevin

run_mcmc_mh_proposal_langevin <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal_langevin.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_langevin",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis-Hastings proposal in MCMC
# ilike::barker

run_mcmc_mh_proposal_barker <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal_barker.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_barker",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in Metropolis-Hastings proposal in MCMC
# ilike::mirror

run_mcmc_mh_proposal_mirror <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_mh_proposal_mirror.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_mh_proposal_mirror",seed=2,initial_values = list(list(theta=1)))
}

# Test custom unadjusted proposal in MCMC

run_mcmc_unadjusted_proposal <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_unadjusted_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_unadjusted_proposal",seed=2,initial_values = list(list(theta=1)))
}

# Test built-in unadjusted proposal in MCMC
# ilike::langevin

run_mcmc_unadjusted_proposal_langevin <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_unadjusted_proposal_langevin.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_unadjusted_proposal_langevin",seed=2,initial_values = list(list(theta=1)))
}

# Test a deterministic sweep: a Metropolis random walk followed by a
# Metropolis-Hastings Langevin move, both used in every iteration.

run_mcmc_sweep_deterministic <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_sweep_deterministic.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_sweep_deterministic",seed=2,initial_values = list(list(theta=1)))
}

# Test a stochastic sweep: the same two moves, one chosen at random in each
# iteration with probabilities given by mcmc_weights.

run_mcmc_sweep_stochastic <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_sweep_stochastic.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_sweep_stochastic",seed=2,initial_values = list(list(theta=1)))
}

# Metropolis-within-Gibbs on a model with two independent parameters: each
# move lists only the factors that depend on the parameter it changes.

run_mcmc_gibbs_m_proposal <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior_two_params.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_gibbs_m_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_gibbs_m_proposal",seed=2,initial_values = list(list(theta1=1,theta2=1)))
}

run_mcmc_gibbs_mh_proposal_langevin <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior_two_params.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_gibbs_mh_proposal_langevin.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_gibbs_mh_proposal_langevin",seed=2,initial_values = list(list(theta1=1,theta2=1)))
}

# Cut model (Plummer 2015): the phi move cuts the likelihood of y.

run_mcmc_cut <- function()
{
  filename_model = system.file("extdata", "test_cut_model.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_cut_mcmc.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_cut",seed=2,initial_values = list(list(phi=1,theta=1)))
}

# Sweep: unadjusted Langevin move, then an MH move using the same factors
# (not allowed: make_mcmc stops).

run_mcmc_sweep_unadjusted_then_mh <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_sweep_unadjusted_then_mh.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_sweep_unadjusted_then_mh",seed=2,initial_values = list(list(theta=1)))
}

# Undeclared cut: run with the consistency check on (should stop), and the
# Gibbs example with the check on (should run).

run_mcmc_undeclared_cut_checked <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_m_proposal_undeclared_cut.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_undeclared_cut",seed=2,initial_values = list(list(theta=1)),check_factor_index=TRUE)
}

run_mcmc_gibbs_m_proposal_checked <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior_two_params.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_gibbs_m_proposal.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_gibbs_m_proposal_checked",seed=2,initial_values = list(list(theta1=1,theta2=1)),check_factor_index=TRUE)
}

# Unadjusted move for theta1 and Metropolis move for theta2, using different
# factors (allowed).

run_mcmc_gibbs_unadjusted_then_m <- function()
{
  filename_model = system.file("extdata", "test_gaussian_posterior_two_params.ilike", package = "ilike")
  filename_mcmc = system.file("extdata", "test_gibbs_unadjusted_then_m.ilike", package = "ilike")
  recipe = ilike::compile(c(filename_model, filename_mcmc))
  ilike::MCMC(recipe,results_name="mcmc_gibbs_unadjusted_then_m",seed=2,initial_values = list(list(theta1=1,theta2=1)))
}
