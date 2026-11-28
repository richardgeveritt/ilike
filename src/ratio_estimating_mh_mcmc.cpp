#include <cmath>
#include "ratio_estimating_mh_mcmc.h"
#include "distributions.h"
#include "particle_tuning.h"
#include "mcmc_adaptor.h"
#include "metropolis_hastings_standard_mcmc_output.h"
#include "proposal_kernel.h"

namespace ilike
{

RatioEstimatingMHMCMC::RatioEstimatingMHMCMC()
  :MetropolisHastingsMCMC()
{
}

RatioEstimatingMHMCMC::RatioEstimatingMHMCMC(size_t number_of_iterations_in,
                                         ProposalKernel* proposal_in,
                                         RatioStep* step_in)
  :MetropolisHastingsMCMC(number_of_iterations_in, proposal_in)
{
  this->ratio = EstimatedRatio(step_in);
}

RatioEstimatingMHMCMC::RatioEstimatingMHMCMC(MCMCTermination* termination_in,
                                         ProposalKernel* proposal_in,
                                         RatioStep* step_in)
  :MetropolisHastingsMCMC(termination_in, proposal_in)
{
  this->ratio = EstimatedRatio(step_in);
}

RatioEstimatingMHMCMC::~RatioEstimatingMHMCMC()
{
}

RatioEstimatingMHMCMC::RatioEstimatingMHMCMC(const RatioEstimatingMHMCMC &another)
  :MetropolisHastingsMCMC(another)
{
  this->make_copy(another);
}

void RatioEstimatingMHMCMC::operator=(const RatioEstimatingMHMCMC &another)
{
  if (this == &another)
    return;
  MetropolisHastingsMCMC::operator=(another);
  this->make_copy(another);
}

Kernel* RatioEstimatingMHMCMC::duplicate() const
{
  return new RatioEstimatingMHMCMC(*this);
}

MCMC* RatioEstimatingMHMCMC::mcmc_duplicate() const
{
  return new RatioEstimatingMHMCMC(*this);
}

MetropolisHastingsMCMC* RatioEstimatingMHMCMC::metropolis_hastings_mcmc_duplicate() const
{
  return new RatioEstimatingMHMCMC(*this);
}

void RatioEstimatingMHMCMC::make_copy(const RatioEstimatingMHMCMC &another)
{
  this->ratio = another.ratio;
}

StandardMCMCOutput* RatioEstimatingMHMCMC::initialise_mcmc_output() const
{
  if (this->adapts())
    return new MetropolisHastingsStandardMCMCOutput(this->metropolis_hastings_mcmc_duplicate());
  return new MetropolisHastingsStandardMCMCOutput(const_cast<RatioEstimatingMHMCMC*>(this), false);
}

bool RatioEstimatingMHMCMC::decide(RandomNumberGenerator &rng,
                          double log_u,
                          const Particle &particle,
                          Particle &proposed_particle,
                          double &acceptance_probability,
                          Particle* refreshed,
                          bool* refreshed_made) const
{
  acceptance_probability = 0.0;
  if (this->ratio.changes_state())
  {
    // (C311) a state step: it makes the new state and gives the target part
    // of the ratio; the proposal's densities added, the reverse one read
    // from the new state
    double target = this->ratio.state_log_term(rng, this->index.get(), particle, proposed_particle, refreshed, refreshed_made);
    if ((target == -arma::datum::inf) || std::isnan(target))
      return false;
    double log_ratio = target + this->proposal->evaluate_kernel(particle, proposed_particle) -
      this->proposal->evaluate_kernel(proposed_particle, particle);
    if (log_ratio >= 0.0)
      acceptance_probability = 1.0;
    else if (!std::isnan(log_ratio))
      acceptance_probability = exp(log_ratio);
    return log_u < log_ratio;
  }
  double numerator = proposed_particle.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
    return false;
  double log_ratio = numerator - denominator +
    this->proposal->evaluate_kernel(particle, proposed_particle) -
    this->proposal->evaluate_kernel(proposed_particle, particle) +
    this->ratio.log_term(rng, particle.parameters, proposed_particle.parameters);
  if (log_ratio >= 0.0)
    acceptance_probability = 1.0;
  else if (!std::isnan(log_ratio))
    acceptance_probability = exp(log_ratio);
  return log_u < log_ratio;
}

void RatioEstimatingMHMCMC::accept(Particle &proposed_particle) const
{
  proposed_particle.set_acceptance(this->proposal.get(),true);
  this->refresh_cut_factors(proposed_particle);
  this->check_carried_over_factors(proposed_particle, this->index.get());
}

Particle RatioEstimatingMHMCMC::move(RandomNumberGenerator &rng,
                            const Particle &particle) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    this->get_new_outputs_index(this->index.get()));
  double log_u = log(runif(rng));
  double acceptance_probability;
  Particle refreshed;
  bool refreshed_made = false;
  bool accepted = this->decide(rng, log_u, particle, proposed_particle, acceptance_probability,
                               this->ratio.refreshes() ? &refreshed : NULL, &refreshed_made);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), acceptance_probability);
  if (accepted)
  {
    this->accept(proposed_particle);
    return proposed_particle;
  }
  // (C312) a rejection with a refreshed current state (MHAAR-RB)
  if (refreshed_made)
    return refreshed;
  return particle;
}

void RatioEstimatingMHMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle) const
{
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    this->get_new_outputs_index(this->index.get()));
  double log_u = log(runif(rng));
  double acceptance_probability;
  Particle refreshed;
  bool refreshed_made = false;
  bool accepted = this->decide(rng, log_u, particle, proposed_particle, acceptance_probability,
                               this->ratio.refreshes() ? &refreshed : NULL, &refreshed_made);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), acceptance_probability);
  if (accepted)
  {
    this->accept(proposed_particle);
    particle = std::move(proposed_particle);
  }
  else if (refreshed_made)
    particle = std::move(refreshed); // (C312)
}

void RatioEstimatingMHMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle,
                                 Particle &scratch) const
{
  this->proposal->move_into(rng,
                            particle,
                            this->get_new_outputs_index(this->index.get()),
                            scratch);
  double log_u = log(runif(rng));
  double acceptance_probability;
  Particle refreshed;
  bool refreshed_made = false;
  bool accepted = this->decide(rng, log_u, particle, scratch, acceptance_probability,
                               this->ratio.refreshes() ? &refreshed : NULL, &refreshed_made);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, scratch, this->proposal.get(), acceptance_probability);
  if (accepted)
  {
    this->accept(scratch);
    particle.swap(scratch);
  }
  else if (refreshed_made)
    particle = std::move(refreshed); // (C312)
}

void RatioEstimatingMHMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                          Particle &particle,
                                          Particle &scratch)
{
  // as move_in_place(rng, particle, scratch), then the adaptation step
  this->proposal->prepare_mcmc_adaptation();
  this->proposal->move_into(rng,
                            particle,
                            this->get_new_outputs_index(this->index.get()),
                            scratch);
  double log_u = log(runif(rng));
  double acceptance_probability;
  Particle refreshed;
  bool refreshed_made = false;
  bool accepted = this->decide(rng, log_u, particle, scratch, acceptance_probability,
                               this->ratio.refreshes() ? &refreshed : NULL, &refreshed_made);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, scratch, this->proposal.get(), acceptance_probability);
  if (accepted)
  {
    this->accept(scratch);
    particle.swap(scratch);
  }
  else if (refreshed_made)
    particle = std::move(refreshed); // (C312)
  // after an acceptance the previous state is in scratch
  AdaptationStep step;
  step.previous = accepted ? &scratch : &particle;
  step.proposed = accepted ? &particle : &scratch;
  step.current = &particle;
  step.acceptance_probability = acceptance_probability;
  this->proposal->mcmc_adapt(step);
}

Particle RatioEstimatingMHMCMC::subsample_move(RandomNumberGenerator &rng,
                                      const Particle &particle) const
{
  Rcpp::stop("Moves with estimated acceptance ratios (exchange) are not available with subsampling.");
}

}
