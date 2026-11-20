#include <cmath>
#include "distributions.h"
#include "utils.h"
#include "metropolis_hastings_mcmc.h"
#include "normal_block_proposal_kernel.h"
#include "particle_tuning.h"
#include "mcmc_adaptor.h"
#include "gaussian_random_walk_proposal_kernel.h"
#include "metropolis_hastings_standard_mcmc_output.h"

namespace ilike
{
MetropolisHastingsMCMC::MetropolisHastingsMCMC()
:MCMC()
{
  this->index.reset();
  this->proposal.reset();
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               ProposalKernel* proposal_in)
:MCMC(number_of_iterations_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               const std::string &variable_name_in,
                                               const arma::mat &proposal_covariance_in)
:MCMC(number_of_iterations_in)
{
  this->index.reset();
  std::vector<std::string> variable_names_in;
  variable_names_in.push_back(variable_name_in);
  std::vector<arma::mat> proposal_covariances_in;
  proposal_covariances_in.push_back(proposal_covariance_in);
  this->proposal.reset(new GaussianRandomWalkProposalKernel(variable_names_in,
                                                        proposal_covariances_in));
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(size_t number_of_iterations_in,
                                               const std::vector<std::string> &variable_names_in,
                                               const std::vector<arma::mat> &proposal_covariances_in)
:MCMC(number_of_iterations_in)
{
  this->index.reset();
  this->proposal.reset(new GaussianRandomWalkProposalKernel(variable_names_in,
                                                        proposal_covariances_in));
  // default to Gaussian random walk
  //this->proposal.reset(ProposalKernel(EvaluateLogMCMCProposalPtr proposal_evaluate_in,
  //                                SimulateMCMCProposalPtr proposal_simulate_in,
  //                                proposal_variances_in));
}

MetropolisHastingsMCMC::MetropolisHastingsMCMC(MCMCTermination* termination_in,
                                               ProposalKernel* proposal_in)
:MCMC(termination_in)
{
  this->proposal.reset(proposal_in);
  this->index.reset();
}

MetropolisHastingsMCMC::~MetropolisHastingsMCMC()
{
  this->proposal.reset();
  
  this->index.reset();
}

//Copy constructor for the MetropolisHastingsMCMC class.
MetropolisHastingsMCMC::MetropolisHastingsMCMC(const MetropolisHastingsMCMC &another)
:MCMC(another)
{
  this->make_copy(another);
}

void MetropolisHastingsMCMC::operator=(const MetropolisHastingsMCMC &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  this->proposal.reset();
  
  this->index.reset();
  
  MCMC::operator=(another);
  this->make_copy(another);
}

Kernel* MetropolisHastingsMCMC::duplicate() const
{
  return( new MetropolisHastingsMCMC(*this));
}

MCMC* MetropolisHastingsMCMC::mcmc_duplicate() const
{
  return( new MetropolisHastingsMCMC(*this));
}

MetropolisHastingsMCMC* MetropolisHastingsMCMC::metropolis_hastings_mcmc_duplicate() const
{
  MetropolisHastingsMCMC* result = new MetropolisHastingsMCMC(*this);
  return result;
}

void MetropolisHastingsMCMC::make_copy(const MetropolisHastingsMCMC &another)
{
  if (another.proposal!=NULL)
    this->proposal.reset(another.proposal->proposal_kernel_duplicate());
  else
    this->proposal.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
}

Particle MetropolisHastingsMCMC::move(RandomNumberGenerator &rng,
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
    if (particle.tuning != NULL)
      observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), 0.0);
    return particle;
  }
  
  double log_ratio = numerator - denominator +
    this->proposal->evaluate_kernel(particle, proposed_particle) -
    this->proposal->evaluate_kernel(proposed_particle, particle);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), acceptance_from_log_ratio(log_ratio));
  if (log_u < log_ratio)
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

void MetropolisHastingsMCMC::move_in_place(RandomNumberGenerator &rng,
                                 Particle &particle) const
{
  // Same steps and random numbers as MetropolisHastingsMCMC::move; a rejection leaves
  // particle as it is instead of returning a copy of it.
  Particle proposed_particle = this->proposal->move(rng,
                                                    particle,
                                                    this->get_new_outputs_index(this->index.get()));
  double log_u = log(runif(rng));
  double numerator = proposed_particle.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  if ( (numerator==-arma::datum::inf) || (denominator==-arma::datum::inf) )
  {
    if (particle.tuning != NULL)
      observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), 0.0);
    return;
  }
  double log_ratio = numerator - denominator +
    this->proposal->evaluate_kernel(particle, proposed_particle) -
    this->proposal->evaluate_kernel(proposed_particle, particle);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), acceptance_from_log_ratio(log_ratio));
  if (log_u < log_ratio)
  {
    proposed_particle.set_acceptance(this->proposal.get(),true);
    this->refresh_cut_factors(proposed_particle);
    this->check_carried_over_factors(proposed_particle, this->index.get());
    particle = std::move(proposed_particle);
  }
}

Particle MetropolisHastingsMCMC::subsample_move(RandomNumberGenerator &rng,
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
  
  if (log_u < numerator - denominator +
      this->proposal->subsample_evaluate_kernel(particle, proposed_particle) -
      this->proposal->subsample_evaluate_kernel(proposed_particle, particle))
  {
    proposed_particle.set_acceptance(this->proposal.get(),true);
    return proposed_particle;
  }
  else
  {
    return particle;
  }
  
}

void MetropolisHastingsMCMC::smc_adapt(SMCOutput* current_state)
{
  this->proposal->smc_adapt(current_state);
}

void MetropolisHastingsMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  this->proposal->ensemble_adapt(current_state);
}

void MetropolisHastingsMCMC::specific_mcmc_adapt(const Particle &current_particle,
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
void MetropolisHastingsMCMC::set_index(Index* index_in)
{
  this->index.reset();
  this->index.reset((index_in!=NULL) ? index_in->duplicate() : NULL);
  if (this->proposal!=NULL)
    this->proposal->set_index(index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisHastingsMCMC::set_index_if_null(Index* index_in)
{
  if ( (this->index==NULL) && (index_in!=NULL) )
    this->index.reset(index_in->duplicate());
  if (this->proposal!=NULL)
    this->proposal->set_index_if_null(index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisHastingsMCMC::set_cut_index(const Index* cut_index_in)
{
  this->store_cut_index(cut_index_in);
  this->update_index_with_cut(this->index.get());
}

void MetropolisHastingsMCMC::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  if (this->proposal!=NULL)
  {
    this->proposal->set_proposal_parameters(proposal_parameters_in);
  }
}

std::vector<const ProposalKernel*> MetropolisHastingsMCMC::get_proposals() const
{
  std::vector<const ProposalKernel*> proposals;
  if (this->proposal!=NULL)
  {
    proposals = this->proposal->get_proposals();
  }
  return proposals;
}

void MetropolisHastingsMCMC::move_in_place(RandomNumberGenerator &rng,
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
  {
    if (particle.tuning != NULL)
      observe_tuned_proposal(particle, scratch, this->proposal.get(), 0.0);
    return;
  }
  double log_ratio = numerator - denominator +
    this->proposal->evaluate_kernel(particle, scratch) -
    this->proposal->evaluate_kernel(scratch, particle);
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, scratch, this->proposal.get(), acceptance_from_log_ratio(log_ratio));
  if (log_u < log_ratio)
  {
    scratch.set_acceptance(this->proposal.get(),true);
    this->refresh_cut_factors(scratch);
    this->check_carried_over_factors(scratch, this->index.get());
    particle.swap(scratch);
  }
}

bool MetropolisHastingsMCMC::adapts() const
{
  return this->proposal->has_mcmc_adaptor();
}

StandardMCMCOutput* MetropolisHastingsMCMC::initialise_mcmc_output() const
{
  // A move that adapts gives each chain its own copy of the kernel; otherwise
  // all the particles' outputs share this one (used read-only by the moves).
  if (this->adapts())
    return new MetropolisHastingsStandardMCMCOutput(this->metropolis_hastings_mcmc_duplicate());
  return new MetropolisHastingsStandardMCMCOutput(const_cast<MetropolisHastingsMCMC*>(this), false);
}

void MetropolisHastingsMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
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
    double log_ratio = numerator - denominator +
      this->proposal->evaluate_kernel(particle, scratch) -
      this->proposal->evaluate_kernel(scratch, particle);
    if (log_ratio >= 0.0)
      acceptance_probability = 1.0;
    else if (!std::isnan(log_ratio))
      acceptance_probability = exp(log_ratio);
    if (particle.tuning != NULL)
      observe_tuned_proposal(particle, scratch, this->proposal.get(), acceptance_probability);
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
  step.rng = &rng;
  this->proposal->mcmc_adapt(step);
}

// (C226) hmc's diagnostics; not for gibbs on a vector, whose kernel is
// built on hmc's but makes no trajectory
static const HMCProposalKernel* diagnosed_hmc(const ProposalKernel* proposal)
{
  const HMCProposalKernel* kernel = dynamic_cast<const HMCProposalKernel*>(proposal);
  if ((kernel == NULL) || (dynamic_cast<const NormalBlockProposalKernel*>(proposal) != NULL))
    return NULL;
  return kernel;
}

void MetropolisHastingsMCMC::diagnostic_names(std::vector<std::string> &names) const
{
  if (diagnosed_hmc(this->proposal.get()) != NULL)
    names = {"leapfrog_steps", "accept_stat", "energy_error", "divergent"};
}

bool MetropolisHastingsMCMC::warming_up() const
{
  return this->proposal->mcmc_warming_up();
}

void MetropolisHastingsMCMC::diagnostic_move_in_place(RandomNumberGenerator &rng,
                                                      Particle &particle,
                                                      Particle &scratch,
                                                      bool adapting,
                                                      double* row)
{
  // as move_in_place(rng, particle, scratch), or adaptive_move_in_place when
  // adapting (the same random numbers), recording one row
  const HMCProposalKernel* kernel = diagnosed_hmc(this->proposal.get());
  if (adapting)
    this->proposal->prepare_mcmc_adaptation();
  double steps = double(kernel->get_steps());
  if (const double* tuning = kernel->tuning_values(particle))
    steps = std::max(1.0, std::round(tuning[1]));
  this->proposal->move_into(rng,
                            particle,
                            this->get_new_outputs_index(this->index.get()),
                            scratch);
  double log_u = log(runif(rng));
  double numerator = scratch.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  bool accepted = false;
  double acceptance_probability = 0.0;
  double log_ratio = -arma::datum::inf;
  if ( (numerator!=-arma::datum::inf) && (denominator!=-arma::datum::inf) )
  {
    log_ratio = numerator - denominator +
      this->proposal->evaluate_kernel(particle, scratch) -
      this->proposal->evaluate_kernel(scratch, particle);
    if (log_ratio >= 0.0)
      acceptance_probability = 1.0;
    else if (!std::isnan(log_ratio))
      acceptance_probability = exp(log_ratio);
    if (particle.tuning != NULL)
      observe_tuned_proposal(particle, scratch, this->proposal.get(), acceptance_probability);
    if (log_u < log_ratio)
    {
      scratch.set_acceptance(this->proposal.get(),true);
      this->refresh_cut_factors(scratch);
      this->check_carried_over_factors(scratch, this->index.get());
      particle.swap(scratch);
      accepted = true;
    }
  }
  else if (particle.tuning != NULL)
    observe_tuned_proposal(particle, scratch, this->proposal.get(), 0.0);
  double energy_error = -log_ratio;
  row[0] = steps;
  row[1] = acceptance_probability;
  row[2] = energy_error;
  row[3] = (energy_error <= 1000.0) ? 0.0 : 1.0; // (above 1000, or not finite)
  if (!adapting)
    return;
  AdaptationStep step;
  step.previous = accepted ? &scratch : &particle;
  step.proposed = accepted ? &particle : &scratch;
  step.current = &particle;
  step.acceptance_probability = acceptance_probability;
  step.rng = &rng;
  this->proposal->mcmc_adapt(step);
}

void MetropolisHastingsMCMC::append_adaptation_state(std::vector<std::string> *names,
                                 std::vector<double> *values,
                                 bool final,
                                 const std::string &prefix) const
{
  this->proposal->append_adaptation_state(names, values, final, prefix);
}

}
