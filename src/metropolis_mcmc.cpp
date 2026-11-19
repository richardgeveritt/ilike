#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "metropolis_mcmc.h"
#include "transform.h"
#include "mcmc_adaptor.h"
#include "symmetric_proposal_kernel.h"
#include "ensemble_kalman_output.h"
#include "metropolis_standard_mcmc_output.h"

namespace ilike
{
MetropolisMCMC::MetropolisMCMC()
:MCMC()
{
  this->index.reset();
  this->proposal.reset();
}

MetropolisMCMC::MetropolisMCMC(size_t number_of_iterations_in,
                               SymmetricProposalKernel* proposal_in)
:MCMC(number_of_iterations_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

MetropolisMCMC::MetropolisMCMC(size_t number_of_iterations_in,
                               const std::vector<Parameters> &initial_points_in,
                               const Parameters &proposal_variances_in)
:MCMC(number_of_iterations_in)
{
  this->index.reset();
  this->proposal.reset();
  // default to Gaussian random walk
  //this->proposal.reset(ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in));
}

MetropolisMCMC::MetropolisMCMC(MCMCTermination* termination_in,
                               SymmetricProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

MetropolisMCMC::~MetropolisMCMC()
{
  this->proposal.reset();
  
  this->index.reset();
}

//Copy constructor for the MetropolisMCMC class.
MetropolisMCMC::MetropolisMCMC(const MetropolisMCMC &another)
:MCMC(another)
{
  this->make_copy(another);
}

void MetropolisMCMC::operator=(const MetropolisMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal.reset();
  
  this->index.reset();
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* MetropolisMCMC::duplicate() const
{
  return( new MetropolisMCMC(*this));
}

MCMC* MetropolisMCMC::mcmc_duplicate() const
{
  return( new MetropolisMCMC(*this));
}

MetropolisMCMC* MetropolisMCMC::metropolis_mcmc_duplicate() const
{
  MetropolisMCMC* result = new MetropolisMCMC(*this);
  return result;
}

void MetropolisMCMC::make_copy(const MetropolisMCMC &another)
{
  if (another.proposal!=NULL)
    this->proposal.reset(another.proposal->symmetric_proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

// the log acceptance ratio: for a proposal on a transformed scale (C209), the
// symmetric step is in z, so the ratio is that of the densities of z,
// pi(x) |dx/dz|, i.e. plus log|dz/dx| at the old point minus at the new
double MetropolisMCMC::log_ratio(double numerator,
                                 double denominator,
                                 const Particle &old_particle,
                                 const Particle &proposed_particle) const
{
  const Transform* transform = this->proposal->get_transform();
  if (transform == NULL)
    return numerator - denominator;
  return numerator - denominator + transform->log_abs_jacobian_determinant(old_particle.parameters) -
    transform->log_abs_jacobian_determinant(proposed_particle.parameters);
}

Particle MetropolisMCMC::move(RandomNumberGenerator &rng,
                              const Particle &particle) const
{
  // Only this move's factors get new outputs in the proposed particle; the
  // others are copied with their cached values. They either cancel in the
  // acceptance ratio or (cut factors) are refreshed after acceptance.
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    this->get_new_outputs_index(this->index.get()));
  
  double log_u = log(runif(rng));
  
  double numerator = proposed_particle.evaluate_likelihoods(this->index.get());
  // The denominator is read from the current particle's per-factor cache over
  // this move's own index. (It used to be particle.target_evaluated, a single
  // number holding whatever index was evaluated last, e.g. the enclosing SMC's
  // index or another move's factor index; see B6 in dev/plan.md.)
  double denominator = particle.get_cached_likelihoods(this->index.get());
  
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
  {
    return particle;
  }
  
  if (log_u < this->log_ratio(numerator, denominator, particle, proposed_particle))
  {
    proposed_particle.set_acceptance(this->proposal.get(),true);
    this->refresh_cut_factors(proposed_particle);
    this->check_carried_over_factors(proposed_particle, this->index.get());
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

void MetropolisMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle) const
{
  // Same steps and random numbers as MetropolisMCMC::move; a rejection leaves
  // particle as it is instead of returning a copy of it.
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    this->get_new_outputs_index(this->index.get()));
  double log_u = log(runif(rng));
  double numerator = proposed_particle.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
    return;
  if (log_u < this->log_ratio(numerator, denominator, particle, proposed_particle))
  {
    proposed_particle.set_acceptance(this->proposal.get(),true);
    this->refresh_cut_factors(proposed_particle);
    this->check_carried_over_factors(proposed_particle, this->index.get());
    particle = std::move(proposed_particle);
  }
}

Particle MetropolisMCMC::subsample_move(RandomNumberGenerator &rng,
                                        const Particle &particle) const
{
  Particle proposed_particle = this->proposal->subsample_move(rng,
                                                              particle);
  
  double log_u = log(runif(rng));
  
  double numerator = proposed_particle.subsample_evaluate_likelihoods(this->index.get());
  double denominator = particle.subsample_target_evaluated;
  
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
  {
    return particle;
  }
  
  if (log_u < this->log_ratio(numerator, denominator, particle, proposed_particle))
  {
    proposed_particle.set_acceptance(this->proposal.get(),true);
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

void MetropolisMCMC::smc_adapt(SMCOutput* current_state)
{
  proposal->smc_adapt(current_state);
}

void MetropolisMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  this->proposal->ensemble_adapt(current_state);
}

void MetropolisMCMC::specific_mcmc_adapt(const Particle &current_particle,
                                         size_t iteration_counter)
{
  this->proposal->mcmc_adapt(current_particle,
                             iteration_counter);
}

// Ownership rule for indices: set_index and set_index_if_null never take
// ownership of index_in. This object keeps its own duplicate (deleted in the
// destructor), and the caller remains responsible for index_in. Previously the
// pointer itself was stored, so when one index was handed to several moves (e.g.
// the moves of a sweep) each destructor deleted it -> double free.
void MetropolisMCMC::set_index(Index* index_in)
{
  this->index.reset();
  this->index.reset((index_in!=NULL) ? index_in->duplicate() : NULL);
  if (this->proposal!=NULL)
    this->proposal->set_index(index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisMCMC::set_index_if_null(Index* index_in)
{
  if ( (this->index==NULL) && (index_in!=NULL) )
    this->index.reset(index_in->duplicate());
  if (this->proposal!=NULL)
    this->proposal->set_index_if_null(index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisMCMC::set_cut_index(const Index* cut_index_in)
{
  this->store_cut_index(cut_index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  if (this->proposal!=NULL)
  {
    this->proposal->set_proposal_parameters(proposal_parameters_in);
  }
}

std::vector<const ProposalKernel*> MetropolisMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals;
  if (this->proposal!=NULL)
  {
    proposals = this->proposal->get_proposals();
  }
  return proposals;
}

void MetropolisMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle,
                                 Particle &scratch) const
{
  // As move_in_place(rng, particle), with the proposal built in scratch.
  this->proposal->move_into(rng,
                            particle,
                            this->get_new_outputs_index(this->index.get()),
                            scratch);
  double log_u = log(runif(rng));
  double numerator = scratch.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
    return;
  if (log_u < this->log_ratio(numerator, denominator, particle, scratch))
  {
    scratch.set_acceptance(this->proposal.get(),true);
    this->refresh_cut_factors(scratch);
    this->check_carried_over_factors(scratch, this->index.get());
    particle.swap(scratch);
  }
}

bool MetropolisMCMC::adapts() const
{
  return this->proposal->has_mcmc_adaptor();
}

StandardMCMCOutput* MetropolisMCMC::initialise_mcmc_output() const
{
  // A move that adapts gives each chain its own copy of the kernel; otherwise
  // all the particles' outputs share this one (used read-only by the moves).
  if (this->adapts())
    return new MetropolisStandardMCMCOutput(this->metropolis_mcmc_duplicate());
  return new MetropolisStandardMCMCOutput(const_cast<MetropolisMCMC*>(this), false);
}

void MetropolisMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                Particle &particle,
                                Particle &scratch)
{
  // As move_in_place(rng, particle, scratch) (same random numbers), and then
  // the adaptation step (dev/plan.md 3.2): alpha = min(1, ratio), 0 when
  // either density is 0.
  this->proposal->prepare_mcmc_adaptation();
  this->proposal->move_into(rng,
                            particle,
                            this->get_new_outputs_index(this->index.get()),
                            scratch);
  double log_u = log(runif(rng));
  double numerator = scratch.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  bool accepted = false;
  double acceptance_probability = 0.0;
  if ( (numerator!=-arma::datum::inf) && (denominator!=-arma::datum::inf) )
  {
    double log_ratio = this->log_ratio(numerator, denominator, particle, scratch);
    if (log_ratio >= 0.0)
      acceptance_probability = 1.0;
    else if (!std::isnan(log_ratio))
      acceptance_probability = exp(log_ratio);
    if (log_u < log_ratio)
    {
      scratch.set_acceptance(this->proposal.get(),true);
      this->refresh_cut_factors(scratch);
      this->check_carried_over_factors(scratch, this->index.get());
      particle.swap(scratch);
      accepted = true;
    }
  }
  // after an acceptance the previous state is in scratch
  AdaptationStep step;
  step.previous = accepted ? &scratch : &particle;
  step.proposed = accepted ? &particle : &scratch;
  step.current = &particle;
  step.acceptance_probability = acceptance_probability;
  this->proposal->mcmc_adapt(step);
}

void MetropolisMCMC::append_adaptation_state(std::vector<std::string> *names,
                                 std::vector<double> *values,
                                 bool final,
                                 const std::string &prefix) const
{
  this->proposal->append_adaptation_state(names, values, final, prefix);
}

}
