#include <cmath>
#include <algorithm>
#include "involutive_mcmc.h"
#include "distributions.h"
#include "particle_tuning.h"
#include "metropolis_hastings_standard_mcmc_output.h"

namespace ilike
{

// ---- InvolutionProposalKernel ------------------------------------------------

InvolutionProposalKernel::InvolutionProposalKernel()
  :ProposalKernel()
{
}

InvolutionProposalKernel::~InvolutionProposalKernel()
{
}

InvolutionProposalKernel::InvolutionProposalKernel(const std::vector<std::string> &state_variables_in,
                                                   InvolutionAuxiliaryPtr auxiliary_in,
                                                   InvolutionLogDensityPtr log_auxiliary_in,
                                                   InvolutionMapPtr map_in,
                                                   InvolutionLogDensityPtr log_jacobian_in,
                                                   size_t check_in)
  :ProposalKernel()
{
  this->state_variables = state_variables_in;
  this->auxiliary = auxiliary_in;
  this->log_auxiliary = log_auxiliary_in;
  this->map = map_in;
  this->log_jacobian = log_jacobian_in;
  this->check = check_in;
}

InvolutionProposalKernel::InvolutionProposalKernel(const InvolutionProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void InvolutionProposalKernel::operator=(const InvolutionProposalKernel &another)
{
  if (this == &another)
    return;
  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* InvolutionProposalKernel::duplicate() const
{
  return new InvolutionProposalKernel(*this);
}

ProposalKernel* InvolutionProposalKernel::proposal_kernel_duplicate() const
{
  return new InvolutionProposalKernel(*this);
}

void InvolutionProposalKernel::make_copy(const InvolutionProposalKernel &another)
{
  this->state_variables = another.state_variables;
  this->auxiliary = another.auxiliary;
  this->log_auxiliary = another.log_auxiliary;
  this->map = another.map;
  this->log_jacobian = another.log_jacobian;
  this->check = another.check;
}

Particle InvolutionProposalKernel::build(const Particle &particle, const Parameters &state, const Index* new_outputs_index) const
{
  return this->proposed_from(particle, state, new_outputs_index);
}

double InvolutionProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                          const Particle &old_particle) const
{
  Rcpp::stop("involution(…) is a move of its own, not a proposal.");
}

double InvolutionProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                    const Particle &old_particle) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

Parameters InvolutionProposalKernel::simulate(RandomNumberGenerator &rng,
                                              const Particle &particle) const
{
  Rcpp::stop("involution(…) is a move of its own, not a proposal.");
}

Parameters InvolutionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                        const Particle &particle) const
{
  return this->simulate(rng, particle);
}

Parameters InvolutionProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                        const std::string &variable,
                                                        const Particle &particle) const
{
  return this->simulate(rng, particle);
}

arma::mat InvolutionProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                             const Particle &proposed_particle,
                                                             const Particle &old_particle)
{
  Rcpp::stop("involution(…): the gradient of a proposal's log density is not available.");
}

arma::mat InvolutionProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                       const Particle &proposed_particle,
                                                                       const Particle &old_particle)
{
  Rcpp::stop("involution(…): the gradient of a proposal's log density is not available.");
}

void InvolutionProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* InvolutionProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> InvolutionProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void InvolutionProposalKernel::set_index(Index* index_in)
{
}

void InvolutionProposalKernel::set_index_if_null(Index* index_in)
{
}

bool InvolutionProposalKernel::can_be_evaluated() const
{
  return true;
}

void InvolutionProposalKernel::set_data(Data* data_in)
{
}

// ---- InvolutiveMCMC -----------------------------------------------------------

InvolutiveMCMC::InvolutiveMCMC()
  :MetropolisHastingsMCMC(), involution(NULL), data(NULL), checked(0)
{
}

InvolutiveMCMC::InvolutiveMCMC(size_t number_of_iterations_in, InvolutionProposalKernel* proposal_in, Data* data_in)
  :MetropolisHastingsMCMC(number_of_iterations_in, proposal_in), involution(proposal_in), data(data_in), checked(0)
{
}

InvolutiveMCMC::InvolutiveMCMC(MCMCTermination* termination_in, InvolutionProposalKernel* proposal_in, Data* data_in)
  :MetropolisHastingsMCMC(termination_in, proposal_in), involution(proposal_in), data(data_in), checked(0)
{
}

InvolutiveMCMC::~InvolutiveMCMC()
{
}

InvolutiveMCMC::InvolutiveMCMC(const InvolutiveMCMC &another)
  :MetropolisHastingsMCMC(another), checked(0)
{
  this->make_copy(another);
}

void InvolutiveMCMC::operator=(const InvolutiveMCMC &another)
{
  if (this == &another)
    return;
  MetropolisHastingsMCMC::operator=(another);
  this->make_copy(another);
}

Kernel* InvolutiveMCMC::duplicate() const
{
  return new InvolutiveMCMC(*this);
}

MCMC* InvolutiveMCMC::mcmc_duplicate() const
{
  return new InvolutiveMCMC(*this);
}

MetropolisHastingsMCMC* InvolutiveMCMC::metropolis_hastings_mcmc_duplicate() const
{
  return new InvolutiveMCMC(*this);
}

void InvolutiveMCMC::make_copy(const InvolutiveMCMC &another)
{
  // (the base class copied the proposal: the involution is that copy)
  this->involution = dynamic_cast<InvolutionProposalKernel*>(this->proposal.get());
  this->data = another.data;
  this->checked = another.checked.load();
}

StandardMCMCOutput* InvolutiveMCMC::initialise_mcmc_output() const
{
  if (this->adapts())
    return new MetropolisHastingsStandardMCMCOutput(this->metropolis_hastings_mcmc_duplicate());
  return new MetropolisHastingsStandardMCMCOutput(const_cast<InvolutiveMCMC*>(this), false);
}

static double largest_difference(const arma::mat &a, const arma::mat &b)
{
  if (arma::size(a) != arma::size(b))
    return arma::datum::inf;
  if (a.n_elem == 0)
    return 0.0;
  return arma::abs(a - b).max() / (1.0 + std::max(arma::abs(a).max(), arma::abs(b).max()));
}

bool InvolutiveMCMC::step(RandomNumberGenerator &rng, const Particle &particle, Particle &proposed_particle) const
{
  const InvolutionProposalKernel &k = *this->involution;
  const Data &data = *this->data;
  // v ~ q(. | x), and xi = (x, v)
  Parameters auxiliary = k.auxiliary(rng, particle.parameters, data);
  std::vector<std::string> auxiliary_names = auxiliary.get_vector_variables();
  Parameters xi = particle.parameters.deep_copy();
  const std::vector<std::string> existing = particle.parameters.get_vector_variables();
  for (const std::string &name : auxiliary_names)
  {
    if (std::find(existing.begin(), existing.end(), name) != existing.end())
      Rcpp::stop("involution(…): the auxiliary variable " + name + " has the name of a variable of the model.");
    xi[name] = auxiliary[name];
  }
  double log_q = k.log_auxiliary(xi, data);
  double log_jacobian = (k.log_jacobian != NULL) ? k.log_jacobian(xi, data) : 0.0;
  // (x', v') = phi(x, v)
  Parameters image = k.map(xi, data);
  Parameters state;
  for (const std::string &name : k.state_variables)
    state[name] = image[name];
  proposed_particle = k.build(particle, state, this->get_new_outputs_index(this->index.get()));
  Parameters xi_image = proposed_particle.parameters.deep_copy();
  for (const std::string &name : auxiliary_names)
    xi_image[name] = image[name];
  // the check that phi is an involution, on the first proposals
  if (this->checked.load() < k.check)
  {
    this->checked++;
    Parameters back = k.map(xi_image, data);
    for (const std::string &name : k.state_variables)
      if (largest_difference(back[name], xi[name]) > 1e-8)
        Rcpp::stop("involution(…): the map is not an involution: applying it twice does not give back " + name + ".");
    for (const std::string &name : auxiliary_names)
      if (largest_difference(back[name], xi[name]) > 1e-8)
        Rcpp::stop("involution(…): the map is not an involution: applying it twice does not give back the auxiliary variable " + name + ".");
  }
  double log_u = log(runif(rng));
  double numerator = proposed_particle.evaluate_likelihoods(this->index.get());
  double denominator = particle.get_cached_likelihoods(this->index.get());
  double acceptance_probability = 0.0;
  bool accepted = false;
  if ((numerator != -arma::datum::inf) && (denominator != -arma::datum::inf))
  {
    double log_ratio = numerator - denominator + k.log_auxiliary(xi_image, data) - log_q + log_jacobian;
    if (log_ratio >= 0.0)
      acceptance_probability = 1.0;
    else if (!std::isnan(log_ratio))
      acceptance_probability = exp(log_ratio);
    accepted = (log_u < log_ratio);
  }
  if (particle.tuning != NULL)
    observe_tuned_proposal(particle, proposed_particle, this->proposal.get(), acceptance_probability);
  if (accepted)
  {
    proposed_particle.set_acceptance(this->proposal.get(), true);
    this->refresh_cut_factors(proposed_particle);
    this->check_carried_over_factors(proposed_particle, this->index.get());
  }
  return accepted;
}

Particle InvolutiveMCMC::move(RandomNumberGenerator &rng, const Particle &particle) const
{
  Particle proposed_particle;
  if (this->step(rng, particle, proposed_particle))
    return proposed_particle;
  return particle;
}

void InvolutiveMCMC::move_in_place(RandomNumberGenerator &rng, Particle &particle) const
{
  Particle proposed_particle;
  if (this->step(rng, particle, proposed_particle))
    particle = std::move(proposed_particle);
}

void InvolutiveMCMC::move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch) const
{
  if (this->step(rng, particle, scratch))
    particle.swap(scratch);
}

void InvolutiveMCMC::adaptive_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch)
{
  Rcpp::stop("involution(…) moves do not adapt.");
}

Particle InvolutiveMCMC::subsample_move(RandomNumberGenerator &rng, const Particle &particle) const
{
  Rcpp::stop("involution(…) moves are not available with subsampling.");
}

}
