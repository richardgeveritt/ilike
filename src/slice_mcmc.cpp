#include "slice_mcmc.h"
#include "slice_point_proposal_kernel.h"
#include "distributions.h"
#include "index.h"
#include "factor_variables.h"
#include "factors.h"

namespace ilike
{

SliceMCMC::SliceMCMC()
:MetropolisHastingsMCMC()
{
}

SliceMCMC::SliceMCMC(size_t number_of_iterations_in, SlicePointProposalKernel* proposal_in)
:MetropolisHastingsMCMC(number_of_iterations_in, proposal_in)
{
}

SliceMCMC::SliceMCMC(MCMCTermination* termination_in, SlicePointProposalKernel* proposal_in)
:MetropolisHastingsMCMC(termination_in, proposal_in)
{
}

SliceMCMC::~SliceMCMC()
{
}

SliceMCMC::SliceMCMC(const SliceMCMC &another)
:MetropolisHastingsMCMC(another)
{
}

void SliceMCMC::operator=(const SliceMCMC &another)
{
  if (this == &another)
    return;
  MetropolisHastingsMCMC::operator=(another);
}

Kernel* SliceMCMC::duplicate() const
{
  return new SliceMCMC(*this);
}

MCMC* SliceMCMC::mcmc_duplicate() const
{
  return new SliceMCMC(*this);
}

MetropolisHastingsMCMC* SliceMCMC::metropolis_hastings_mcmc_duplicate() const
{
  return new SliceMCMC(*this);
}

Particle SliceMCMC::move(RandomNumberGenerator &rng,
                         const Particle &particle) const
{
  const SlicePointProposalKernel* point = static_cast<const SlicePointProposalKernel*>(this->proposal.get());
  const Index* new_outputs_index = this->get_new_outputs_index(this->index.get());
  const double lower = point->get_lower(), upper = point->get_upper(), width = point->get_width();
  const size_t max_steps = point->get_max_steps();
  const bool whole_interval = !(width > 0.0);

  // (the slice is a level of the exact target: an estimated one would not
  // leave the target invariant)
  if ((particle.factor_variables != NULL) && (particle.factor_variables->get_factors() != NULL) &&
      particle.factor_variables->get_factors()->any_estimate())
    Rcpp::stop("slice: slice sampling needs the exact target, but a factor's value is an estimate (an inner "
               "algorithm, a synthetic likelihood or an estimated normalising constant). Use a Metropolis or "
               "Metropolis-Hastings move.");
  double log_target = particle.get_cached_likelihoods(this->index.get());
  if (log_target == -arma::datum::inf)
    return particle; // (no slice: the current point has density 0)
  Particle current = particle;
  arma::colvec x = arma::vectorise(arma::mat(current.parameters[point->get_variable()]));
  bool moved = false;

  // the log target at x with element k set to v (-Inf outside [lower, upper]),
  // and the particle there
  auto evaluate = [&](arma::uword k, double v, Particle &candidate) -> double
  {
    if ((v < lower) || (v > upper))
      return -arma::datum::inf;
    point->set_point(k, v);
    candidate = point->move(rng, current, new_outputs_index);
    return candidate.evaluate_likelihoods(this->index.get());
  };

  Particle candidate;
  for (arma::uword k = 0; k < x.n_elem; ++k)
  {
    double log_y = log_target + std::log(runif(rng));
    double x0 = x[k];
    double L, R;
    if (whole_interval)
    {
      L = lower;
      R = upper;
    }
    else
    {
      // stepping out (Neal 2003, fig. 3), within [lower, upper]
      L = x0 - width*runif(rng);
      R = L + width;
      size_t J = size_t(std::floor(double(max_steps)*runif(rng)));
      size_t K = (max_steps - 1) - std::min(J, max_steps - 1);
      while ((J > 0) && (L > lower) && (evaluate(k, L, candidate) > log_y))
      {
        L = L - width;
        --J;
      }
      while ((K > 0) && (R < upper) && (evaluate(k, R, candidate) > log_y))
      {
        R = R + width;
        --K;
      }
      L = std::max(L, lower);
      R = std::min(R, upper);
    }
    // shrinkage (fig. 5): x0 is on the slice, so the interval cannot shrink
    // to nothing; the cap guards against a target that is not exact
    for (size_t s = 0; s < 1000; ++s)
    {
      double x1 = L + runif(rng)*(R - L);
      double log_target_1 = evaluate(k, x1, candidate);
      if (log_target_1 > log_y)
      {
        current = std::move(candidate);
        log_target = log_target_1;
        x[k] = x1;
        moved = true;
        break;
      }
      if (x1 < x0)
        L = x1;
      else
        R = x1;
    }
  }
  if (!moved)
    return particle;
  current.set_acceptance(this->proposal.get(), true);
  this->refresh_cut_factors(current);
  return current;
}

void SliceMCMC::move_in_place(RandomNumberGenerator &rng,
                              Particle &particle) const
{
  particle = this->move(rng, particle);
}

void SliceMCMC::move_in_place(RandomNumberGenerator &rng,
                              Particle &particle,
                              Particle &scratch) const
{
  particle = this->move(rng, particle);
}

void SliceMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                       Particle &particle,
                                       Particle &scratch)
{
  particle = this->move(rng, particle);
}

Particle SliceMCMC::subsample_move(RandomNumberGenerator &rng,
                                   const Particle &particle) const
{
  Rcpp::stop("slice with subsampling is not available.");
}

}
