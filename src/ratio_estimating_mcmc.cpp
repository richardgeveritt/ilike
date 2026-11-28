#include <cmath>
#include "ratio_estimating_mcmc.h"
#include "distributions.h"
#include "particle_tuning.h"
#include "mcmc_adaptor.h"
#include "symmetric_proposal_kernel.h"
#include "metropolis_standard_mcmc_output.h"

namespace ilike
{

RatioEstimatingMCMC::RatioEstimatingMCMC()
  :MetropolisMCMC()
{
}

RatioEstimatingMCMC::RatioEstimatingMCMC(size_t number_of_iterations_in,
                                         SymmetricProposalKernel* proposal_in,
                                         RatioStep* step_in)
  :MetropolisMCMC(number_of_iterations_in, proposal_in)
{
  this->ratio = EstimatedRatio(step_in);
}

RatioEstimatingMCMC::RatioEstimatingMCMC(MCMCTermination* termination_in,
                                         SymmetricProposalKernel* proposal_in,
                                         RatioStep* step_in)
  :MetropolisMCMC(termination_in, proposal_in)
{
  this->ratio = EstimatedRatio(step_in);
}

RatioEstimatingMCMC::~RatioEstimatingMCMC()
{
}

RatioEstimatingMCMC::RatioEstimatingMCMC(const RatioEstimatingMCMC &another)
  :MetropolisMCMC(another)
{
  this->make_copy(another);
}

void RatioEstimatingMCMC::operator=(const RatioEstimatingMCMC &another)
{
  if (this == &another)
    return;
  MetropolisMCMC::operator=(another);
  this->make_copy(another);
}

Kernel* RatioEstimatingMCMC::duplicate() const
{
  return new RatioEstimatingMCMC(*this);
}

MCMC* RatioEstimatingMCMC::mcmc_duplicate() const
{
  return new RatioEstimatingMCMC(*this);
}

MetropolisMCMC* RatioEstimatingMCMC::metropolis_mcmc_duplicate() const
{
  return new RatioEstimatingMCMC(*this);
}

void RatioEstimatingMCMC::make_copy(const RatioEstimatingMCMC &another)
{
  this->ratio = another.ratio;
}

StandardMCMCOutput* RatioEstimatingMCMC::initialise_mcmc_output() const
{
  if (this->adapts())
    return new MetropolisStandardMCMCOutput(this->metropolis_mcmc_duplicate());
  return new MetropolisStandardMCMCOutput(const_cast<RatioEstimatingMCMC*>(this), false);
}

bool RatioEstimatingMCMC::decide(RandomNumberGenerator &rng,
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
    // of the ratio; the proposal's own terms (a transform's Jacobian) added
    double target = this->ratio.state_log_term(rng, this->index.get(), particle, proposed_particle, refreshed, refreshed_made);
    if ((target == -arma::datum::inf) || std::isnan(target))
      return false;
    double log_ratio = this->log_ratio(0.0, 0.0, particle, proposed_particle) + target;
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
  double log_ratio = this->log_ratio(numerator, denominator, particle, proposed_particle) +
    this->ratio.log_term(rng, particle.parameters, proposed_particle.parameters);
  if (log_ratio >= 0.0)
    acceptance_probability = 1.0;
  else if (!std::isnan(log_ratio))
    acceptance_probability = exp(log_ratio);
  return log_u < log_ratio;
}

void RatioEstimatingMCMC::accept(Particle &proposed_particle) const
{
  proposed_particle.set_acceptance(this->proposal.get(),true);
  this->refresh_cut_factors(proposed_particle);
  this->check_carried_over_factors(proposed_particle, this->index.get());
}

Particle RatioEstimatingMCMC::move(RandomNumberGenerator &rng,
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

void RatioEstimatingMCMC::move_in_place(RandomNumberGenerator &rng,
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

void RatioEstimatingMCMC::move_in_place(RandomNumberGenerator &rng,
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

void RatioEstimatingMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
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

Particle RatioEstimatingMCMC::subsample_move(RandomNumberGenerator &rng,
                                      const Particle &particle) const
{
  Rcpp::stop("Moves with estimated acceptance ratios (exchange) are not available with subsampling.");
}

}
