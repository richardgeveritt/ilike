// The table of built-ins (dev/plan.md 3.14, I2): one entry per built-in,
// with its arguments and its constructors. Adding a built-in means adding an
// entry here (and writing the class it builds). The binding of arguments,
// keyword arguments, run-time constants and the messages are in
// builtin_registry.cpp.

#include <RcppArmadillo.h>
#include <cctype>
#include <cmath>

#include "builtin_registry.h"

#include "gaussian_distribution_factor.h"
#include "loggaussian_distribution_factor.h"
#include "gamma_distribution_factor.h"
#include "gaussian_independent_proposal_kernel.h"
#include "loggaussian_independent_proposal_kernel.h"
#include "gamma_independent_proposal_kernel.h"
#include "gaussian_random_walk_proposal_kernel.h"
#include "uniform_random_walk_proposal_kernel.h"
#include "componentwise_random_walk_proposal_kernel.h"
#include "principal_components_proposal_kernel.h"
#include "langevin_proposal_kernel.h"
#include "hmc_proposal_kernel.h"
#include "normal_block_proposal_kernel.h"
#include "conjugate_block_proposal_kernel.h"
#include "dual_averaging_adaptor.h"
#include "pretune_adaptor.h"
#include "truncated_normal_distribution_factor.h"
#include "truncated_normal_independent_proposal_kernel.h"
#include "truncated_normal_random_walk_proposal_kernel.h"
#include "dimension_jump_proposal_kernel.h"
#include "pcn_proposal_kernel.h"
#include "pcn_langevin_proposal_kernel.h"
#include "slice_point_proposal_kernel.h"
#include "adaptive_iterations_mcmc_termination.h"
#include "esjd_adaptor.h"
#include "barker_dynamics_proposal_kernel.h"
#include "mirror_proposal_kernel.h"
#include "direct_gradient_estimator.h"
#include "gaussian_mcmc_adaptor.h"
#include "proposal_parameter_adaptor.h"
#include "gaussian_smc_adaptor.h"
#include "iterations_mcmc_termination.h"
#include "ess_smc_criterion.h"
#include "cess_smc_criterion.h"
#include "positive_smc_criterion.h"
#include "stable_smc_termination.h"
#include "always_smc_termination.h"
#include "stochastic_ensemble_shifter.h"
#include "square_root_ensemble_shifter.h"
#include "adjustment_ensemble_shifter.h"

namespace ilike
{

typedef ArgumentType T;
typedef BuiltinRole R;

static const std::vector<BuiltinRole> distribution_roles = {R::Prior, R::ISProposal, R::IMHProposal};
static const std::vector<BuiltinRole> random_walk_roles = {R::MHProposal, R::UnadjustedProposal, R::MProposal,
                                                           R::TransitionModel, R::TransitionProposal};

static ArgumentSpec the_variable()
{
  return required("variable", T::Variable, "the variable");
}

// --- MCMC adaptors (dev/plan.md 3.2; Andrieu & Thoms 2008) --------------------

// The default target acceptance rate for a proposal (0: by the variable's
// dimension, for the mirror proposal: mirror_default_target_acceptance).
static double default_target(const std::string &proposal_type)
{
  std::string p = proposal_type;
  for (char &c : p)
    c = std::tolower((unsigned char)c);
  if (p == "componentwise_rw" || p == "pc_rw")
    return 0.44;
  if (p == "mirror")
    return 0.0;
  if (p == "langevin")
    return 0.574;
  if (p == "barker")
    return 0.4;
  if ((p == "hmc") || (p == "nuts"))
    return 0.8; // (owner's decision, dev/plan.md 3.18)
  return 0.234;
}

static bool is_one_of(const std::string &type, const std::vector<std::string> &names)
{
  for (const std::string &name : names)
  {
    std::string a = type, b = name;
    if (a.size() != b.size())
      continue;
    bool same = true;
    for (size_t i = 0; i < a.size(); ++i)
      if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i]))
        same = false;
    if (same)
      return true;
  }
  return false;
}

// which proposals an adaptor applies to, for the check and its message
struct AdaptorFit
{
  std::vector<std::string> proposals;
  bool user_written;
  std::string fits;
};

static void check_fit(const BuiltinArguments &arguments, const MCMCAdaptorContext &context,
                      const AdaptorFit &fit)
{
  const std::string &proposal_type = context.proposal_type;
  bool user_written = proposal_type.empty();
  bool ok = user_written ? fit.user_written : is_one_of(proposal_type, fit.proposals);
  if (ok)
    return;
  std::string label = context.block + ",ilike::" + arguments.builtin();
  if (is_one_of(proposal_type, {"componentwise_rw"}))
    Rcpp::stop(label + ": the proposal ilike::componentwise_rw is adapted by ilike::componentwise_am.");
  if (is_one_of(proposal_type, {"pc_rw"}))
    Rcpp::stop(label + ": the proposal ilike::pc_rw is adapted by ilike::pc_am or ilike::online_pc_am.");
  Rcpp::stop(label + ": applies to " + fit.fits + ".");
}

// target (if the method has one), kappa, C; the gain C (n+1)^(-kappa)
static void adaptor_settings(const BuiltinArguments &arguments, const MCMCAdaptorContext &context,
                             double &target, double &kappa, double &constant)
{
  target = default_target(context.proposal_type);
  if (arguments.given("target"))
    target = arguments.real("target");
  kappa = arguments.real("kappa");
  constant = arguments.real("C");
  if (constant * std::pow(2.0, -kappa) >= 1.0)
    Rcpp::stop(context.block + ": the gain C (n+1)^(-kappa) must be below 1 from the "
               "first update (C 2^(-kappa) < 1).");
}

static std::function<MCMCAdaptor*(const BuiltinArguments&, const MCMCAdaptorContext&)>
gaussian_adaptor(GaussianAdaptationMethod method, AdaptorFit fit)
{
  return [method, fit](const BuiltinArguments &arguments, const MCMCAdaptorContext &context) -> MCMCAdaptor* {
    check_fit(arguments, context, fit);
    double target, kappa, constant;
    adaptor_settings(arguments, context, target, kappa, constant);
    size_t refresh = arguments.has("refresh") ? size_t(arguments.integer("refresh")) : 0;
    std::vector<std::string> variables(1, arguments.variable("variable"));
    return new GaussianMCMCAdaptor(variables, method, target, kappa, constant, refresh,
                                   arguments.logical("averaged"));
  };
}

static ArgumentSpec adaptor_target()
{
  return optional("target", T::Real,
                  "the target acceptance rate (default: 0.234 for random walks and named parameters, "
                  "0.44 for componentwise and principal-components moves, by dimension for the mirror "
                  "proposal, 0.574 for Langevin, 0.4 for Barker)").between(0.0, 1.0);
}
static ArgumentSpec adaptor_kappa()
{
  return optional("kappa", T::Real, RecipeValue(0.6), "the exponent of the gain C (n+1)^(-kappa)");
}
static ArgumentSpec adaptor_constant()
{
  return optional("C", T::Real, RecipeValue(1.0), "the constant of the gain C (n+1)^(-kappa)");
}
static ArgumentSpec adaptor_refresh()
{
  return optional("refresh", T::Integer, RecipeValue(0.0),
                  "update the proposal every refresh iterations (0: the default, every d iterations; "
                  "every iteration for online_pc_am)").at_least(0.0);
}
static ArgumentSpec adaptor_averaged()
{
  return flag("averaged", "use the averages of the iterates (Polyak & Juditsky 1992)");
}

static const AdaptorFit am_fit = {{"norm_rw", "mvnorm_rw", "mirror", "hmc", "nuts"}, false,
                                  "the proposals ilike::norm_rw, ilike::mvnorm_rw, ilike::mirror, ilike::hmc and ilike::nuts"};
static const AdaptorFit scale_fit = {{"norm_rw", "mvnorm_rw", "mirror", "langevin", "barker", "hmc", "nuts"}, false,
                                     "the proposals ilike::norm_rw, ilike::mvnorm_rw, ilike::mirror, "
                                     "ilike::langevin, ilike::barker, ilike::hmc and ilike::nuts"};
static const AdaptorFit componentwise_fit = {{"componentwise_rw"}, false, "the proposal ilike::componentwise_rw"};
static const AdaptorFit pc_fit = {{"pc_rw"}, false, "the proposal ilike::pc_rw"};

// --- the table ---------------------------------------------------------------

static std::vector<Builtin> make_registry()
{
  std::vector<Builtin> table;

  // distributions: priors, and independent proposals
  {
    Builtin b;
    b.name = "norm";
    b.description = "normal distribution";
    b.arguments = {the_variable(), required("mean", T::Real, "the mean"),
                   required("sd", T::Real, "the standard deviation")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new GaussianDistributionFactor(a.variable("variable"), a.real("mean"), a.real("sd")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new GaussianIndependentProposalKernel(a.variable("variable"), a.real("mean"), a.real("sd")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "mvnorm";
    b.description = "multivariate normal distribution";
    b.arguments = {the_variable(), required("mean", T::Vector, "the mean"),
                   required("cov", T::Matrix, "the covariance matrix")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new GaussianDistributionFactor(a.variable("variable"), a.vector("mean"), a.matrix("cov")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new GaussianIndependentProposalKernel(a.variable("variable"), a.vector("mean"), a.matrix("cov")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "lnorm";
    b.description = "log-normal distribution";
    b.arguments = {the_variable(), required("meanlog", T::Real, "the mean of the logarithm"),
                   required("sdlog", T::Real, "the standard deviation of the logarithm")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new LogGaussianDistributionFactor(a.variable("variable"), a.real("meanlog"), a.real("sdlog")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new LogGaussianIndependentProposalKernel(a.variable("variable"), a.real("meanlog"), a.real("sdlog")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "mvlnorm";
    b.description = "multivariate log-normal distribution";
    b.arguments = {the_variable(), required("meanlog", T::Vector, "the mean of the logarithm"),
                   required("covlog", T::Matrix, "the covariance matrix of the logarithm")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new LogGaussianDistributionFactor(a.variable("variable"), a.vector("meanlog"), a.matrix("covlog")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new LogGaussianIndependentProposalKernel(a.variable("variable"), a.vector("meanlog"), a.matrix("covlog")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "gamma";
    b.description = "gamma distribution";
    b.arguments = {the_variable(), required("shape", T::Real, "the shape"),
                   required("rate", T::Real, "the rate")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new GammaDistributionFactor(a.variable("variable"), a.real("shape"), a.real("rate")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new GammaIndependentProposalKernel(a.variable("variable"), a.real("shape"), a.real("rate")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "tnorm";
    b.description = "truncated normal distribution (C217)";
    b.arguments = {the_variable(), required("mean", T::Real, "the mean of the normal"),
                   required("sd", T::Real, "the standard deviation of the normal").above(0.0),
                   optional("lower", T::Real, RecipeValue(-arma::datum::inf), "the lower limit (default -Inf)"),
                   optional("upper", T::Real, RecipeValue(arma::datum::inf), "the upper limit (default Inf)")};
    b.roles = distribution_roles;
    b.prior = [](const BuiltinArguments &a) -> DistributionFactor* {
      return new TruncatedNormalDistributionFactor(a.variable("variable"), a.real("mean"), a.real("sd"),
                                                   a.real("lower"), a.real("upper")); };
    b.independent_proposal = [](const BuiltinArguments &a) -> IndependentProposalKernel* {
      return new TruncatedNormalIndependentProposalKernel(a.variable("variable"), a.real("mean"), a.real("sd"),
                                                          a.real("lower"), a.real("upper")); };
    table.push_back(b);
  }

  // random walks
  {
    Builtin b;
    b.name = "tnorm_rw";
    b.description = "truncated normal random walk: each element to a normal(x, sd) draw truncated to [lower, upper] (C217)";
    b.arguments = {the_variable(), required("sd", T::Real, "the standard deviation of the step").above(0.0),
                   optional("lower", T::Real, RecipeValue(-arma::datum::inf), "the lower limit (default -Inf)"),
                   optional("upper", T::Real, RecipeValue(arma::datum::inf), "the upper limit (default Inf)")};
    // (not symmetric: the truncated mass depends on the current point)
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new TruncatedNormalRandomWalkProposalKernel(a.variable("variable"), a.real("sd"), a.real("lower"), a.real("upper")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "dimension_jump";
    b.description = "a model indicator to k + 1 or k - 1 (forced at the bounds): the proposal of a reversible-jump move (birth_death; C313)";
    b.arguments = {the_variable(), required("minimum", T::Real, "the smallest value of the indicator"),
                   required("maximum", T::Real, "the largest value of the indicator")};
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new DimensionJumpProposalKernel(a.variable("variable"), a.real("minimum"), a.real("maximum")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "pcn";
    b.description = "preconditioned Crank-Nicolson: u* = mean + rho (u - mean) + sqrt(1 - rho^2) N(0, cov); "
                    "leaves N(mean, cov) invariant (dev/plan.md 3.21)";
    b.arguments = {the_variable(), required("mean", T::Vector, "the mean of u's Gaussian distribution"),
                   required("cov", T::Matrix, "the covariance of u's Gaussian distribution"),
                   required("rho", T::Real, "the autocorrelation, in [0, 1)").at_least(0.0)};
    // (not symmetric: q(u* | u) is N(mean + rho (u - mean), (1 - rho^2) cov))
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new PCNProposalKernel(a.variable("variable"), a.vector("mean"), a.matrix("cov"), a.real("rho")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "pcn_mala";
    b.description = "preconditioned Crank-Nicolson Langevin (pCN-MALA, Cotter et al. 2013): pCN for u with distribution "
                    "N(mean, cov), with a step along the gradient of the other factors; delta > 0 the step size";
    b.arguments = {the_variable(), required("mean", T::Vector, "the mean of u's Gaussian distribution"),
                   required("cov", T::Matrix, "the covariance of u's Gaussian distribution"),
                   required("delta", T::Real, "the step size (delta > 0; small: close to u)").above(0.0)};
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new PCNLangevinProposalKernel(a.variable("variable"), a.vector("mean"), a.matrix("cov"), a.real("delta"),
                                           new DirectGradientEstimator()); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "slice";
    b.description = "slice sampling, each element in turn (a move: move m = slice(...); Neal 2003): the interval "
                    "is [lower, upper] when no width is given, otherwise stepped out by width";
    b.arguments = {the_variable(),
                   optional("width", T::Real, RecipeValue(0.0), "the step of the stepping out (0: the whole of [lower, upper])").at_least(0.0),
                   optional("lower", T::Real, RecipeValue(-arma::datum::inf), "the lower limit (default -Inf)"),
                   optional("upper", T::Real, RecipeValue(arma::datum::inf), "the upper limit (default Inf)"),
                   optional("max_steps", T::Integer, RecipeValue(50.0), "the most steps of the stepping out").at_least(1.0)};
    // (the settings of a SliceMCMC made by the recipe, as nuts's)
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new SlicePointProposalKernel(a.variable("variable"), a.real("width"), a.real("lower"), a.real("upper"),
                                          size_t(a.integer("max_steps"))); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "norm_rw";
    b.description = "normal random walk";
    b.arguments = {the_variable(), required("sd", T::Real, "the standard deviation of the step")};
    b.roles = random_walk_roles;
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new GaussianRandomWalkProposalKernel(a.variable("variable"), a.real("sd")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "mvnorm_rw";
    b.description = "multivariate normal random walk, step N(0, scale cov)";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the covariance of the step"),
                   optional("scale", T::Real, RecipeValue(1.0), "the factor s^2 multiplying the covariance")};
    b.roles = random_walk_roles;
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new GaussianRandomWalkProposalKernel(a.variable("variable"), a.matrix("cov"), a.real("scale")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "unif_rw";
    b.description = "uniform random walk, step U(-halfwidth, halfwidth)";
    b.arguments = {the_variable(), required("halfwidth", T::Real, "the half-width of the step")};
    b.roles = random_walk_roles;
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new UniformRandomWalkProposalKernel(a.variable("variable"), a.real("halfwidth")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "mvunif_rw";
    b.description = "multivariate uniform random walk";
    b.arguments = {the_variable(), required("halfwidth", T::Vector, "the half-widths of the step")};
    b.roles = random_walk_roles;
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      arma::colvec halfwidth = a.vector("halfwidth");
      return new UniformRandomWalkProposalKernel(a.variable("variable"), halfwidth); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "componentwise_rw";
    b.description = "random-scan componentwise normal random walk, component k moved by N(0, scale cov[k,k]) "
                    "(Andrieu & Thoms 2008, Algorithm 5)";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the covariance (its diagonal is used)"),
                   optional("scale", T::Real, RecipeValue(1.0), "the factor s^2 multiplying the variances")};
    b.roles = {R::MProposal};
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new ComponentwiseRandomWalkProposalKernel(a.variable("variable"), a.matrix("cov"), a.real("scale")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "pc_rw";
    b.description = "random walk along one principal direction of cov at a time "
                    "(Andrieu & Thoms 2008, Algorithm 8)";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the covariance"),
                   optional("scale", T::Real, RecipeValue(1.0), "the factor s^2 multiplying the covariance"),
                   optional("directions", T::Integer, RecipeValue(0.0),
                            "use the first directions principal directions (0: all)").at_least(0.0)};
    b.roles = {R::MProposal};
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new PrincipalComponentsProposalKernel(a.variable("variable"), a.matrix("cov"), a.real("scale"),
                                                   size_t(a.integer("directions"))); };
    table.push_back(b);
  }

  // gradient-based and other proposals
  {
    Builtin b;
    b.name = "langevin";
    b.description = "Langevin (MALA) proposal";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the preconditioning covariance"),
                   required("scale", T::Real, "the factor s^2 multiplying the covariance")};
    b.roles = {R::MHProposal, R::UnadjustedProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new LangevinProposalKernel(a.variable("variable"), a.matrix("cov"), a.real("scale"),
                                        new DirectGradientEstimator()); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "hmc";
    b.description = "Hamiltonian Monte Carlo proposal: steps leapfrog steps of size step, momentum N(0, mass^-1)";
    b.arguments = {the_variable(), required("step", T::Real, "the leapfrog step size epsilon").above(0.0),
                   optional("steps", T::Integer, "the number of leapfrog steps L").at_least(1.0),
                   optional("length", T::Real, "the trajectory length (L = max(1, round(length/step)))").above(0.0),
                   optional("mass", T::Matrix, "the inverse mass matrix (a covariance; default: the identity)")};
    b.roles = {R::MHProposal, R::UnadjustedProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      if (a.given("steps") == a.given("length"))
        Rcpp::stop("ilike::hmc: give the number of leapfrog steps (steps = L) or the trajectory length (length = …), one of them.");
      double step = a.real("step");
      size_t steps = a.given("steps") ? size_t(a.integer("steps"))
                                      : size_t(std::max(1.0, std::round(a.real("length")/step)));
      arma::mat mass = a.given("mass") ? a.matrix("mass") : arma::mat(1, 1, arma::fill::ones);
      return new HMCProposalKernel(a.variable("variable"), mass, step, steps); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "nuts";
    b.description = "the No-U-Turn sampler (a move: move m = nuts(...)), leapfrog step size step, momentum N(0, mass^-1)";
    b.arguments = {the_variable(), required("step", T::Real, "the leapfrog step size epsilon").above(0.0),
                   optional("mass", T::Matrix, "the inverse mass matrix (a covariance; default: the identity)"),
                   optional("max_depth", T::Integer, RecipeValue(10.0), "the maximum depth of the tree (at most 2^max_depth - 1 leapfrog steps)").at_least(1.0),
                   optional("variant", T::String, RecipeValue(std::string("stan")),
                            "stan (the default: multinomial sampling, the generalised U-turn criterion) or paper (Hoffman & Gelman's Algorithm 3: slice sampling)")};
    // (the role of its Hamiltonian system, in a NoUTurnMCMC made by the recipe)
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      arma::mat mass = a.given("mass") ? a.matrix("mass") : arma::mat(1, 1, arma::fill::ones);
      HMCProposalKernel* kernel = new HMCProposalKernel(a.variable("variable"), mass, a.real("step"), 1);
      kernel->set_max_depth(size_t(a.integer("max_depth")));
      std::string variant = a.string("variant");
      if ((variant != "paper") && (variant != "stan"))
        Rcpp::stop("ilike::nuts: variant is paper or stan.");
      kernel->set_stan_variant(variant == "stan");
      return kernel; };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "normal_block";
    b.description = "gibbs on a vector (a move: move m = gibbs(beta)): a draw from the block's normal full conditional";
    b.arguments = {the_variable()};
    // (the proposal of the Metropolis-Hastings move that gibbs on a vector is
    // lowered to, C216; refused when written directly)
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new NormalBlockProposalKernel(a.variable("variable")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "dirichlet_block";
    b.description = "gibbs on a Dirichlet-distributed probability vector (a move: move m = gibbs(p))";
    b.arguments = {the_variable()};
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new DirichletBlockProposalKernel(a.variable("variable")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "normal_invgamma_block";
    b.description = "gibbs on a mean and a variance (a move: move m = gibbs([mu, v]))";
    b.arguments = {the_variable(), optional("mean_first", T::Logical, RecipeValue(true), "whether the mean comes first")};
    b.roles = {R::MHProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new NormalInverseGammaBlockProposalKernel(a.variable("variable"), a.logical("mean_first")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "barker";
    b.description = "Barker proposal";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the preconditioning covariance"),
                   required("scale", T::Real, "the factor s^2 multiplying the covariance")};
    b.roles = {R::MHProposal, R::UnadjustedProposal};
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      return new BarkerDynamicsProposalKernel(a.variable("variable"), a.matrix("cov"), a.real("scale"),
                                              new DirectGradientEstimator()); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "hamiltonian";
    b.description = "Hamiltonian Monte Carlo proposal (not written yet)";
    b.arguments = {the_variable(), required("cov", T::Matrix, "the mass matrix")};
    b.roles = {R::MHProposal, R::UnadjustedProposal};
    // HMCProposalKernel is a stub (simulate returns no parameters, so the
    // chain never moves); refuse it rather than return a stuck chain.
    b.proposal = [](const BuiltinArguments &a) -> ProposalKernel* {
      Rcpp::stop("MH proposal ilike::hamiltonian (HMC) is not written yet.");
      return nullptr; };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "mirror";
    b.description = "mirror proposal: reflection of the current state through mean, plus N(0, scale cov) noise";
    b.arguments = {the_variable(), required("mean", T::Vector, "the centre of the reflection"),
                   required("cov", T::Matrix, "the covariance of the noise"),
                   required("scale", T::Real, "the factor s^2 multiplying the covariance")};
    b.roles = {R::MHProposal, R::UnadjustedProposal, R::MProposal};
    b.symmetric_proposal = [](const BuiltinArguments &a) -> SymmetricProposalKernel* {
      return new MirrorProposalKernel(a.variable("variable"), a.vector("mean"), a.matrix("cov"), a.real("scale")); };
    table.push_back(b);
  }

  // MCMC adaptors
  {
    Builtin b;
    b.name = "scaled_am";
    b.description = "adaptive Metropolis with a scale adapted to a target acceptance rate "
                    "(Andrieu & Thoms 2008, Algorithm 4)";
    b.arguments = {the_variable(), adaptor_target(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(),
                   adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::SCALED_AM, am_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "am";
    b.description = "adaptive Metropolis (Andrieu & Thoms 2008, Algorithm 2)";
    b.arguments = {the_variable(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(), adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::AM, am_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "rb_am";
    b.description = "Rao-Blackwellised adaptive Metropolis (Andrieu & Thoms 2008, Algorithm 3)";
    b.arguments = {the_variable(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(), adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::RB_AM, am_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "scale";
    b.description = "the scale recursion of Andrieu & Thoms (2008), Algorithm 4, alone";
    b.arguments = {the_variable(), adaptor_target(), adaptor_kappa(), adaptor_constant(), adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::SCALE, scale_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "dual_averaging";
    b.description = "the step size of hmc or nuts by dual averaging (Hoffman & Gelman 2014, Algorithms 4-6); needs until = …";
    b.arguments = {the_variable(),
                   optional("target", T::Real, "the target acceptance statistic delta (default 0.8)").between(0.0, 1.0),
                   optional("gamma", T::Real, RecipeValue(0.05), "the shrinkage towards mu = log(10 eps0)").above(0.0),
                   optional("t0", T::Real, RecipeValue(10.0), "stabilises the first iterations").at_least(0.0),
                   optional("kappa", T::Real, RecipeValue(0.75), "the weights of the average of log eps, m^-kappa").between(0.0, 1.0),
                   optional("find_step", T::Logical, RecipeValue(true),
                            "start from FindReasonableEpsilon (Algorithm 4) rather than the step given"),
                   optional("mass", T::String, RecipeValue(std::string("none")),
                            "none, diagonal or dense: also estimate the inverse mass matrix in Stan's windows "
                            "(needs until = iterations(n))")};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = [](const BuiltinArguments &arguments, const MCMCAdaptorContext &context) -> MCMCAdaptor* {
      static const AdaptorFit fit = {{"hmc", "nuts"}, false, "the proposals ilike::hmc and ilike::nuts"};
      check_fit(arguments, context, fit);
      double target = arguments.given("target") ? arguments.real("target") : 0.8;
      std::string mass = arguments.string("mass");
      if ((mass != "none") && (mass != "diagonal") && (mass != "dense"))
        Rcpp::stop(context.block + ": dual_averaging's mass is none, diagonal or dense.");
      return new DualAveragingAdaptor(target, arguments.real("gamma"), arguments.real("t0"), arguments.real("kappa"),
                                      arguments.logical("find_step"), mass); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "componentwise_am";
    b.description = "componentwise adaptive scaling (Andrieu & Thoms 2008, Algorithm 5)";
    b.arguments = {the_variable(), adaptor_target(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(),
                   adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::COMPONENTWISE_AM, componentwise_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "pc_am";
    b.description = "principal-components adaptation by eigen-decomposition every refresh iterations "
                    "(Andrieu & Thoms 2008, Algorithm 8)";
    b.arguments = {the_variable(), adaptor_target(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(),
                   adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::PC_EIGEN, pc_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "online_pc_am";
    b.description = "principal-components adaptation by the online PCA recursion "
                    "(Andrieu & Thoms 2008, Algorithm 8)";
    b.arguments = {the_variable(), adaptor_target(), adaptor_kappa(), adaptor_constant(), adaptor_refresh(),
                   adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = gaussian_adaptor(GaussianAdaptationMethod::PC_ONLINE, pc_fit);
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "acceptance_rate";
    b.description = "adapts a named proposal parameter of a proposal written by the user to a target "
                    "acceptance rate";
    b.arguments = {required("parameter", T::Variable, "the name of the proposal parameter"), adaptor_target(),
                   adaptor_kappa(), adaptor_constant(), adaptor_averaged()};
    b.roles = {R::MCMCAdaptor};
    b.mcmc_adaptor = [](const BuiltinArguments &a, const MCMCAdaptorContext &context) -> MCMCAdaptor* {
      check_fit(a, context, AdaptorFit{{}, true, "proposals written by the user (with proposal_parameters)"});
      double target, kappa, constant;
      adaptor_settings(a, context, target, kappa, constant);
      std::vector<std::string> names(1, a.variable("parameter"));
      return new ProposalParameterAdaptor(names, target, kappa, constant, a.logical("averaged")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "averaged";
    b.description = "Polyak-Juditsky averaging for a user's adaptation function (adapt_m_proposal, "
                    "adapt_mh_proposal)";
    b.roles = {R::MCMCAdaptor};
    // handled where the adaptor is attached (attach_mcmc_adaptor), since it
    // goes with a user's function, not with a built-in adaptor
    b.mcmc_adaptor = [](const BuiltinArguments &a, const MCMCAdaptorContext &context) -> MCMCAdaptor* {
      std::string prefix = context.block.substr(0, context.block.find("_adaptor"));
      Rcpp::stop(context.block + ",ilike::averaged() goes with an adapt_" + prefix + "_proposal function; "
                 "for built-in adaptors, give \"averaged\" as an argument (e.g. ilike::scaled_am(theta,averaged)).");
      return nullptr; };
    table.push_back(b);
  }

  // SMC adaptor
  {
    Builtin b;
    b.name = "esjd";
    b.description = "tunes the proposal for each particle from its squared jump distance (Fearnhead & Taylor 2013; "
                    "Buchholz, Chopin & Jacob 2021): scale (norm_rw, mvnorm_rw, langevin), step and steps (hmc), step (nuts)";
    b.arguments = {the_variable(),
                   optional("scale", T::Vector, "the range of the first multipliers of the step (random walks, langevin)"),
                   optional("step", T::Vector, "the range of the first step sizes (hmc, nuts)"),
                   optional("steps", T::Vector, "the range of the first numbers of leapfrog steps (hmc)"),
                   optional("a", T::Real, RecipeValue(0.0), "the constant a of the weights a + jump").at_least(0.0),
                   optional("perturb", T::Real, RecipeValue(0.05), "the sd of the multiplicative perturbation (log scale)").above(0.0),
                   optional("metric", T::String, RecipeValue(std::string("diagonal")), "diagonal or full: the particles' variances or covariance")};
    b.roles = {R::SMCAdaptor};
    b.smc_adaptor = [](const BuiltinArguments &a, ProposalKernel *proposal) -> SMCAdaptor* {
      std::string variable = a.variable("variable");
      auto range = [&](const std::string &name) {
        arma::colvec r = a.vector(name);
        if ((r.n_elem != 2) || !(r[0] > 0.0) || !(r[0] <= r[1]))
          Rcpp::stop("esjd: " + name + " = [lower, upper], with 0 < lower <= upper.");
        return r;
      };
      std::vector<std::string> names;
      std::vector<double> lower, upper;
      std::vector<bool> integer;
      int cost = -1;
      HMCProposalKernel* hamiltonian = dynamic_cast<HMCProposalKernel*>(proposal);
      if (hamiltonian != NULL)
      {
        if (a.given("scale"))
          Rcpp::stop("esjd: hmc and nuts are tuned by step (and steps for hmc), not scale.");
        if (!a.given("step"))
          Rcpp::stop("esjd: give the range of the step sizes, step = [lower, upper].");
        arma::colvec r = range("step");
        names = {"step", "steps"};
        lower = {r[0], 0.0};
        upper = {r[1], 0.0};
        integer = {false, true};
        if (a.given("steps"))
        {
          arma::colvec s = range("steps");
          lower[1] = std::max(1.0, std::round(s[0]));
          upper[1] = std::max(lower[1], std::round(s[1]));
        }
        else
          lower[1] = upper[1] = double(hamiltonian->get_steps());
        cost = 1; // (a leapfrog step's cost: Buchholz et al.)
      }
      else if ((dynamic_cast<GaussianRandomWalkProposalKernel*>(proposal) != NULL) ||
               (dynamic_cast<LangevinProposalKernel*>(proposal) != NULL))
      {
        if (a.given("step") || a.given("steps"))
          Rcpp::stop("esjd: random walks and langevin are tuned by scale (a multiplier of the step).");
        if (!a.given("scale"))
          Rcpp::stop("esjd: give the range of the first multipliers of the step, scale = [lower, upper].");
        arma::colvec r = range("scale");
        names = {"scale"};
        lower = {r[0]};
        upper = {r[1]};
        integer = {false};
      }
      else
        Rcpp::stop("esjd is available for ilike::norm_rw, ilike::mvnorm_rw, ilike::langevin, ilike::hmc and ilike::nuts.");
      std::string metric = a.string("metric");
      if ((metric != "diagonal") && (metric != "full"))
        Rcpp::stop("esjd: metric is \"diagonal\" or \"full\".");
      return new ESJDAdaptor(variable, names, lower, upper, integer, a.real("a"), a.real("perturb"), metric == "full", cost); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "pretune";
    b.description = "tunes hmc's step and steps for each particle from trial moves before each SMC move "
                    "(Buchholz, Chopin & Jacob 2021, pre-tuning)";
    b.arguments = {the_variable(),
                   optional("max_step", T::Real, "the first upper limit eps* of the trial steps (default: twice hmc's step)").above(0.0),
                   optional("max_steps", T::Integer, RecipeValue(100.0), "the first upper limit L_max of the trial numbers of steps").at_least(1.0),
                   optional("target", T::Real, RecipeValue(0.9), "the acceptance rate that sets eps*").between(0.0, 1.0),
                   optional("metric", T::String, RecipeValue(std::string("diagonal")), "diagonal or full: the particles' variances or covariance"),
                   optional("raise_share", T::Real, RecipeValue(0.5), "L_max is raised when more than this share of the chosen L ...").between(0.0, 1.0),
                   optional("raise_above", T::Real, RecipeValue(0.8), "... are above this fraction of L_max").between(0.0, 1.0),
                   optional("lower_share", T::Real, RecipeValue(0.9), "L_max is lowered when more than this share of the chosen L ...").between(0.0, 1.0),
                   optional("lower_below", T::Real, RecipeValue(0.5), "... are below this fraction of L_max").between(0.0, 1.0),
                   optional("steps_change", T::Integer, RecipeValue(5.0), "how much L_max is raised or lowered (0: fixed)").at_least(0.0)};
    b.roles = {R::SMCAdaptor};
    b.smc_adaptor = [](const BuiltinArguments &a, ProposalKernel *proposal) -> SMCAdaptor* {
      HMCProposalKernel* hamiltonian = dynamic_cast<HMCProposalKernel*>(proposal);
      if (hamiltonian == NULL)
        Rcpp::stop("pretune is available for ilike::hmc (in metropolis_hastings).");
      std::string metric = a.string("metric");
      if ((metric != "diagonal") && (metric != "full"))
        Rcpp::stop("pretune: metric is \"diagonal\" or \"full\".");
      double max_step = a.given("max_step") ? a.real("max_step") : 2.0*hamiltonian->step_size();
      return new PretuneAdaptor(a.variable("variable"), max_step, size_t(a.integer("max_steps")), a.real("target"),
                                metric == "full", a.real("raise_share"), a.real("raise_above"), a.real("lower_share"),
                                a.real("lower_below"), size_t(a.integer("steps_change"))); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "gaussian";
    b.description = "sets the proposal's covariance (and a mirror's centre, a principal-components "
                    "proposal's directions) from the weighted particles before each move";
    b.arguments = {the_variable(),
                   optional("c", T::Real, "with c, the scale is set to c^2/d (c^2 for componentwise and "
                                          "principal-components proposals)").above(0.0)};
    b.roles = {R::SMCAdaptor};
    b.smc_adaptor = [](const BuiltinArguments &a, ProposalKernel *proposal) -> SMCAdaptor* {
      std::string variable = a.variable("variable");
      double constant = a.given("c") ? a.real("c") : -1.0;
      if (proposal->get_gaussian_proposal_info(variable) == NULL)
        Rcpp::stop("smc_adaptor,ilike::gaussian: the proposal has no Gaussian covariance "
                   "for " + variable + " (available for ilike::norm_rw, ilike::mvnorm_rw, "
                   "ilike::componentwise_rw, ilike::pc_rw, ilike::mirror, ilike::langevin, "
                   "ilike::barker, and the independent proposals ilike::norm and ilike::mvnorm).");
      return new GaussianSMCAdaptor(variable, constant); };
    table.push_back(b);
  }

  // MCMC terminations
  {
    Builtin b;
    b.name = "iterations";
    b.description = "a fixed number of iterations";
    b.arguments = {required("n", T::Integer, "the number of iterations").at_least(0.0)};
    b.roles = {R::MCMCTermination};
    b.mcmc_termination = [](const BuiltinArguments &a) -> MCMCTermination* {
      return new IterationsMCMCTermination(size_t(a.integer("n"))); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "adaptive_iterations";
    b.description = "in SMC, each step's number of iterations from the particles' acceptance rate p in one "
                    "pilot iteration: ceil(log c / log(1 - p)), so each particle moves with probability about "
                    "1 - c (South, Pettitt & Drovandi 2019)";
    b.arguments = {optional("c", T::Real, RecipeValue(0.2), "the probability of a particle not moving"),
                   optional("max", T::Integer, RecipeValue(100.0), "the largest number of iterations").at_least(1.0)};
    b.roles = {R::MCMCTermination};
    b.mcmc_termination = [](const BuiltinArguments &a) -> MCMCTermination* {
      return new AdaptiveIterationsMCMCTermination(a.real("c"), size_t(a.integer("max"))); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "se";
    b.description = "stop when the Monte Carlo standard error is below a threshold (not written yet)";
    b.arguments = {required("threshold", T::Real, "the threshold"),
                   optional("max_iterations", T::Integer, "the maximum number of iterations").at_least(0.0)};
    b.roles = {R::MCMCTermination};
    // SEMCMCTermination::terminate() ignores its threshold and stops only at
    // the maximum number of iterations, so refuse it until the rule is written.
    b.mcmc_termination = [](const BuiltinArguments &a) -> MCMCTermination* {
      Rcpp::stop("MCMC termination ilike::se is not written yet: use ilike::iterations.");
      return nullptr; };
    table.push_back(b);
  }

  // SMC: adaptive resampling and adaptive targets (criteria)
  {
    Builtin b;
    b.name = "ess";
    b.description = "effective sample size: resample when the ESS falls below proportion times the number "
                    "of particles (adaptive_resampling), or choose the next target so that the ESS is that "
                    "value (adaptive_target)";
    b.arguments = {required("proportion", T::Real, "the proportion of the number of particles")};
    b.roles = {R::AdaptiveResampling, R::AdaptiveTarget};
    b.smc_criterion = [](const BuiltinArguments &a, size_t number_of_particles) -> SMCCriterion* {
      return new ESSSMCCriterion(a.real("proportion") * double(number_of_particles)); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "cess";
    b.description = "conditional effective sample size: choose the next target so that the CESS is "
                    "proportion times the number of particles";
    b.arguments = {required("proportion", T::Real, "the proportion of the number of particles")};
    b.roles = {R::AdaptiveTarget};
    b.smc_criterion = [](const BuiltinArguments &a, size_t number_of_particles) -> SMCCriterion* {
      return new CESSSMCCriterion(a.real("proportion") * double(number_of_particles)); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "positive";
    b.description = "no adaptation: the next target is the next value of the schedule (the default when none is given)";
    b.roles = {R::AdaptiveTarget};
    b.smc_criterion = [](const BuiltinArguments &a, size_t number_of_particles) -> SMCCriterion* {
      return new PositiveSMCCriterion(); };
    table.push_back(b);
  }

  // SMC terminations
  {
    Builtin b;
    b.name = "stable";
    b.description = "stop when the absolute value of the algorithm's termination score is at most "
                    "threshold for consecutive iterations in a row";
    b.arguments = {required("threshold", T::Real, "the threshold"),
                   required("consecutive", T::Integer, "the number of iterations in a row").at_least(1.0)};
    b.roles = {R::SMCTermination};
    b.smc_termination = [](const BuiltinArguments &a) -> SMCTermination* {
      // the constructor takes (number in a row, threshold)
      return new StableSMCTermination(size_t(a.integer("consecutive")), a.real("threshold")); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "always";
    b.description = "stop at the first check";
    b.roles = {R::SMCTermination};
    b.smc_termination = [](const BuiltinArguments &a) -> SMCTermination* {
      return new AlwaysSMCTermination(); };
    table.push_back(b);
  }

  // ensemble Kalman shifters
  {
    Builtin b;
    b.name = "stochastic";
    b.description = "stochastic (perturbed observations) ensemble shift";
    b.roles = {R::EnKShifter};
    b.enk_shifter = [](const BuiltinArguments &a) -> EnsembleShifter* { return new StochasticEnsembleShifter(); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "sqrt";
    b.description = "square-root ensemble shift";
    b.roles = {R::EnKShifter};
    b.enk_shifter = [](const BuiltinArguments &a) -> EnsembleShifter* { return new SquareRootEnsembleShifter(); };
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "adjustment";
    b.description = "adjustment ensemble shift";
    b.roles = {R::EnKShifter};
    b.enk_shifter = [](const BuiltinArguments &a) -> EnsembleShifter* { return new AdjustmentEnsembleShifter(); };
    table.push_back(b);
  }

  // SMC sequences (the schedules are values: evaluated by the R binding,
  // which runs R expressions such as seq(0, 1, length.out = 11))
  {
    Builtin b;
    b.name = "annealing";
    b.description = "anneal the factors with a power variable along a schedule from 0 to 1";
    b.arguments = {required("variable", T::Variable, "the name of the power"),
                   required("schedule", T::Vector, "the powers"),
                   optional("factors", T::Vector, "the factors (numbered from 1) raised to the power (default: all)")};
    b.roles = {R::SMCSequence, R::EnKSequence};
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "tempering";
    b.description = "as annealing";
    b.arguments = {required("variable", T::Variable, "the name of the power"),
                   required("schedule", T::Vector, "the powers"),
                   optional("factors", T::Vector, "the factors (numbered from 1) raised to the power (default: all)")};
    b.roles = {R::SMCSequence};
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "parameter";
    b.description = "a sequence of values of a parameter of the target (e.g. an ABC tolerance)";
    b.arguments = {required("variable", T::Variable, "the name of the parameter"),
                   required("schedule", T::Vector, "the values"),
                   optional("factors", T::Vector, "the factors (numbered from 1) that depend on the parameter "
                                                  "(default: all)")};
    b.roles = {R::SMCSequence};
    table.push_back(b);
  }

  {
    Builtin b;
    b.name = "data";
    b.description = "data tempering: the data points of factors written as model statements enter along a schedule of "
                    "fractions from 0 to 1, each annealed in";
    b.arguments = {required("variable", T::Variable, "the name of the fraction of the data"),
                   required("schedule", T::Vector, "the fractions"),
                   optional("factors", T::Vector, "the factors (numbered from 1) whose data points enter (default: all "
                                                  "the likelihood factors)"),
                   optional("order", T::String, RecipeValue(std::string("random")),
                            "the order in which the data points enter: random (fixed by the seed) or given"),
                   optional("anneal", T::Logical, RecipeValue(true),
                            "anneal each data point in (else add whole data points)")};
    b.roles = {R::SMCSequence};
    table.push_back(b);
  }

  // simulation-based likelihoods (bound only: built in recipe_factors.cpp)
  {
    auto points = [] { return required("number_of_points", T::Integer, "the number of simulations per evaluation").at_least(0.0); };
    auto tolerance_variable = [] { return required("tolerance_variable", T::Variable,
                                                   "the name of the tolerance (an SMC sequence on it makes it adaptive)"); };
    auto store = [] { return optional("store_output", T::Logical, RecipeValue(false), "store the simulations"); };
    auto scale_variable = [] { return optional("scale_variable", T::Variable, "the name of a variable holding the scale of the summaries"); };
    auto scale = [] { return optional("scale", T::Vector, "the scale of the summaries"); };
    auto parallel = [] { return optional("parallel", T::Logical, RecipeValue(false), "simulate in parallel"); };
    auto grain = [](double g) { return optional("grain_size", T::Integer, RecipeValue(g), "simulations per parallel task").at_least(0.0); };

    Builtin b;
    b.roles = {R::SBILikelihood};
    b.name = "euclidean_uniform_abc";
    b.description = "ABC likelihood with a uniform kernel on the Euclidean distance between summaries";
    b.arguments = {points(), tolerance_variable(), required("tolerance", T::Real, "the tolerance (if not adaptive)"),
                   store(), scale_variable(), scale(), parallel(), grain(100000)};
    table.push_back(b);

    b.name = "gaussian_abc";
    b.description = "ABC likelihood with a Gaussian kernel";
    table.push_back(b);

    b.name = "lp_uniform_abc";
    b.description = "ABC likelihood with a uniform kernel on the L^p distance between summaries";
    b.arguments = {points(), tolerance_variable(), required("tolerance", T::Real, "the tolerance (if not adaptive)"),
                   required("p", T::Real, "the power of the L^p distance"),
                   store(), scale_variable(), scale(), parallel(), grain(100000)};
    table.push_back(b);

    b.name = "ienki_abc";
    b.description = "ABC by iterative ensemble Kalman inversion (epsilon = 0 with estimator_type 0-2: "
                    "synthetic likelihood)";
    b.arguments = {points(), required("estimator_type", T::Integer, "the estimator (0-2: see the documentation)").at_least(0.0),
                   tolerance_variable(), required("epsilon", T::Real, "the tolerance"),
                   required("shifter", T::String, "the ensemble shifter: stochastic, sqrt or adjustment"),
                   required("number_of_targets", T::Integer, "the number of intermediate targets").at_least(1.0),
                   optional("significance_level", T::Real, RecipeValue(1.0), "the significance level"),
                   optional("lag", T::Integer, RecipeValue(0.0), "the lag").at_least(0.0),
                   scale_variable(), scale(),
                   optional("on_summary", T::Logical, RecipeValue(true), "work on the summaries"),
                   parallel(), grain(100000)};
    table.push_back(b);

    b.name = "sl";
    b.description = "synthetic likelihood";
    b.arguments = {points(), optional("unbiased", T::Logical, RecipeValue(false), "use the unbiased estimator"),
                   store(), parallel(), grain(100000)};
    table.push_back(b);
  }

  // algorithmic likelihoods: a filter run on an inner model file (bound only)
  {
    auto model = [] { return required("model", T::String, "the inner model file"); };
    auto parallel = [] { return optional("parallel", T::Logical, RecipeValue(false), "run in parallel"); };
    Builtin b;
    b.roles = {R::AlgorithmicLikelihood};
    b.name = "kf";
    b.description = "Kalman filter on a linear-Gaussian inner model";
    b.arguments = {model(), optional("iterations_to_store", T::Integer, RecipeValue(0.0), "filter iterations stored").at_least(0.0),
                   // (C276) off by default
                   optional("steady_state", T::Real, RecipeValue(0.0),
                            "freeze the covariances once the predicted covariance changes by at most this (relative) between steps; 0: never").at_least(0.0)};
    table.push_back(b);

    b.name = "enkf";
    b.description = "ensemble Kalman filter on an inner model";
    b.arguments = {model(), required("ensemble_members", T::Integer, "the number of ensemble members").at_least(0.0),
                   required("shifter", T::String, "the ensemble shifter: stochastic, sqrt or adjustment"),
                   optional("iterations_to_store", T::Integer, RecipeValue(0.0), "filter iterations stored").at_least(0.0),
                   parallel(),
                   optional("grain_size", T::Integer, RecipeValue(100000.0), "members per parallel task").at_least(0.0),
                   // (dev/plan.md 3.21 (c2'''), C267) set by increase = exchange(below, multiply, max)
                   optional("exchange_below", T::Real, RecipeValue(),
                            "SMC^2: exchange for more members when the acceptance rate is below this").at_least(0.0),
                   optional("exchange_multiply", T::Real, RecipeValue(2.0), "SMC^2: the factor by which an exchange multiplies the members").at_least(1.0),
                   optional("exchange_max", T::Integer, RecipeValue(), "SMC^2: the largest number of members after an exchange").at_least(1.0),
                   // (dev/plan.md 3.22 U9 R3b) correlated pseudo-marginal MCMC
                   optional("correlated", T::Logical, RecipeValue(false),
                            "MCMC: each member's random numbers from its own seed, and a proposal keeping all but refresh members' (correlated estimates)"),
                   optional("refresh", T::Integer, RecipeValue(1.0), "MCMC with correlated = true: the number of members whose random numbers a proposal replaces").at_least(1.0),
                   // (C306) set by correlated = crank_nicolson(sd = s)
                   optional("crank_nicolson_sd", T::Real, RecipeValue(), "MCMC: the Crank-Nicolson step on the members' standard normals").at_least(0.0)};
    table.push_back(b);

    b.name = "pf";
    b.description = "particle filter on an inner model";
    b.arguments = {model(), required("particles", T::Integer, "the number of particles").at_least(0.0),
                   optional("iterations_to_store", T::Integer, RecipeValue(1.0),
                            "generations stored for each outer point (1: the final one)").at_least(0.0),
                   parallel(),
                   optional("grain_size", T::Integer, RecipeValue(64.0), "particles per parallel task").at_least(0.0),
                   // (dev/plan.md 3.21 (c2''), F2) set by increase = exchange(below, multiply, max)
                   optional("exchange_below", T::Real, RecipeValue(),
                            "SMC^2: exchange for more particles when the acceptance rate is below this").at_least(0.0),
                   optional("exchange_multiply", T::Real, RecipeValue(2.0), "SMC^2: the factor by which an exchange multiplies the particles").at_least(1.0),
                   optional("exchange_max", T::Integer, RecipeValue(), "SMC^2: the largest number of particles after an exchange").at_least(1.0),
                   // (dev/plan.md 3.22 U9 R3c, C307) set by correlated = crank_nicolson(sd = s)
                   optional("correlated", T::Logical, RecipeValue(false), "MCMC: correlated estimates (correlated = crank_nicolson(sd = s))"),
                   optional("crank_nicolson_sd", T::Real, RecipeValue(), "MCMC: the Crank-Nicolson step on the filter's standard normals").at_least(0.0)};
    table.push_back(b);

    // integrate u with is(points = N) (dev/plan.md 3.21, step (b)): the
    // "model" is the variable integrated out (set by the lowering)
    b.name = "is";
    b.description = "importance sampling over a variable integrated out of the model";
    b.arguments = {required("model", T::String, "the variable integrated out"),
                   required("points", T::Integer, "the number of importance points").at_least(1.0),
                   parallel(),
                   optional("grain_size", T::Integer, RecipeValue(100000.0), "points per parallel task").at_least(0.0)};
    table.push_back(b);

    // integrate u with smc(particles = N, move = …, sequence = …) (dev/plan.md
    // 3.21, (c1)): the moves, sequence, termination, target and resampling
    // are lowered into the inner algorithm, not passed here
    b.name = "smc";
    b.description = "SMC with MCMC moves over a variable integrated out of the model";
    b.arguments = {required("model", T::String, "the variable integrated out"),
                   required("particles", T::Integer, "the number of particles").at_least(1.0),
                   parallel(),
                   optional("grain_size", T::Integer, RecipeValue(64.0), "particles per parallel task").at_least(0.0)};
    table.push_back(b);
  }

  // state-space model settings (bound only)
  {
    Builtin b;
    b.name = "ssm";
    b.description = "the time index of a state-space model";
    b.arguments = {required("index", T::Variable, "the name of the time index"),
                   required("first_index", T::Integer, "the first index").at_least(0.0),
                   required("last_index", T::Integer, "the last index").at_least(0.0),
                   required("time", T::Variable, "the name of the time variable"),
                   required("initial_time", T::Real, "the initial time"),
                   required("time_difference", T::Variable, "the name of the time-step variable"),
                   required("time_step", T::Real, "the time between updates"),
                   required("predictions_per_update", T::Integer, "the number of predictions per update").at_least(0.0)};
    b.roles = {R::SSM};
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "linear_gaussian_data_model";
    b.description = "data y ~ N(A x, covariance) given the state x";
    b.arguments = {required("variable", T::Variable, "the data variable"),
                   required("state_variables", T::String, "the state variable(s), separated by ;"),
                   required("matrix", T::Matrix, "the matrix A"),
                   required("covariance", T::Matrix, "the covariance")};
    b.roles = {R::LinearGaussianDataModel};
    table.push_back(b);
  }
  {
    Builtin b;
    b.name = "linear_gaussian_transition_model";
    b.description = "state x_t ~ N(A x_{t-1}, covariance)";
    b.arguments = {required("variable", T::Variable, "the state variable"),
                   required("matrix", T::Matrix, "the matrix A"),
                   required("covariance", T::Matrix, "the covariance")};
    b.roles = {R::LinearGaussianTransitionModel};
    table.push_back(b);
  }

  // every role listed has a constructor
  for (const Builtin &b : table)
    for (BuiltinRole role : b.roles)
    {
      bool ok = true;
      switch (role)
      {
        case R::Prior: ok = bool(b.prior) && bool(b.independent_proposal); break;
        case R::ISProposal: case R::IMHProposal: ok = bool(b.independent_proposal); break;
        case R::MHProposal: case R::UnadjustedProposal: ok = bool(b.proposal) || bool(b.symmetric_proposal); break;
        case R::MProposal: case R::TransitionModel: case R::TransitionProposal: ok = bool(b.symmetric_proposal); break;
        case R::MCMCAdaptor: ok = bool(b.mcmc_adaptor); break;
        case R::SMCAdaptor: ok = bool(b.smc_adaptor); break;
        case R::MCMCTermination: ok = bool(b.mcmc_termination); break;
        case R::AdaptiveResampling: case R::AdaptiveTarget: ok = bool(b.smc_criterion); break;
        case R::SMCTermination: ok = bool(b.smc_termination); break;
        case R::EnKShifter: ok = bool(b.enk_shifter); break;
        case R::SMCSequence: case R::EnKSequence: case R::SBILikelihood: case R::AlgorithmicLikelihood: case R::SSM:
        case R::LinearGaussianDataModel: case R::LinearGaussianTransitionModel: break;
      }
      if (!ok)
        Rcpp::stop("Registry: ilike::" + b.name + " lists " + role_block(role) + " but has no constructor for it.");
    }

  return table;
}

const std::vector<Builtin> &builtin_registry()
{
  static const std::vector<Builtin> table = make_registry();
  return table;
}

}
