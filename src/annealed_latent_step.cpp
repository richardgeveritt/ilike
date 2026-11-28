#include <cmath>
#include "annealed_latent_step.h"
#include "particle.h"
#include "index.h"
#include "factor_variables.h"
#include "mcmc.h"

namespace ilike
{

AnnealedLatentStep::AnnealedLatentStep(size_t bridges_in, SMCPath* path_in, Index* state_index_in,
                                       const std::vector<size_t> &kernel_in)
  : annealing(bridges_in, path_in, state_index_in, kernel_in)
{
  if ((bridges_in > 0) && (path_in == NULL))
    Rcpp::stop("annealed(…): the path's SMC was not built.");
}

AnnealedLatentStep::AnnealedLatentStep(const AnnealedLatentStep &another)
  : annealing(another.annealing)
{
}

AnnealedLatentStep::~AnnealedLatentStep()
{
}

RatioStep* AnnealedLatentStep::duplicate() const
{
  return new AnnealedLatentStep(*this);
}

double AnnealedLatentStep::draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const
{
  Rcpp::stop("AnnealedLatentStep::draw: an annealed move's step changes the state (draw_state).");
}

void AnnealedLatentStep::link_kernel(const std::vector<const MCMC*> &cycle)
{
  this->annealing.link_kernel(cycle, "annealed");
}

double AnnealedLatentStep::draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const
{
  // the variables the proposal changed: theta and theta'
  std::vector<std::string> names;
  Parameters start, end;
  {
    const Parameters &before = current.parameters;
    const Parameters &after = proposed.parameters;
    for (auto a = after.vector_begin(); a != after.vector_end(); ++a)
    {
      const arma::mat* b = before.find_vector(a->first);
      if ((b == NULL) || (b->n_elem != a->second.first->n_elem))
        continue;
      if ((b == a->second.first.get()) || arma::approx_equal(*b, *a->second.first, "absdiff", 0.0))
        continue;
      names.push_back(a->first);
      start[a->first] = *b;
      end[a->first] = *a->second.first;
    }
  }
  // (the new state is made on a copy of the current particle)
  Particle annealed = current;
  if (names.empty())
  {
    // (the proposal changed nothing: the new state is the current one)
    proposed = std::move(annealed);
    return 0.0;
  }
  const double cached = annealed.get_cached_likelihoods(index);
  if (cached == -arma::datum::inf)
    return -arma::datum::inf;
  double log_ratio = 0.0;
  if (this->annealing.bridges == 0)
  {
    for (const std::string &name : names)
      annealed.parameters.overwrite_variable(name, end[name]);
  }
  else
  {
    // (A2) not with an estimated factor among those that read theta (U8 replaces this)
    this->annealing.refuse_estimates(current, index, "annealed");
    // the path: the general engine with one particle from the current state
    log_ratio = this->annealing.path->single_state(rng, annealed.parameters, names, start, end, current.tuning);
    if (std::isnan(log_ratio) || (log_ratio == -arma::datum::inf))
      return -arma::datum::inf;
    // (the outputs carried over from the current particle's: the state changed)
    if ((annealed.factor_variables != NULL) && (current.factor_variables != NULL))
      annealed.factor_variables->redo_carry_over(current.factor_variables.get(), annealed.parameters, this->annealing.state_index.get());
  }
  // (cached transforms and gradients belong to the old state: made again for
  // the new one, as a new particle's are)
  annealed.current_proposal_store.clear();
  annealed.simulate_proposal_variables(annealed.proposals_to_transform_for_pointer, annealed.proposals_to_find_gradient_for_pointer);
  // the move's factors at the new state (with no bridges, the whole ratio);
  // after a path, which changed the latent variables too, every factor (the
  // particle's cached values, read by the next moves)
  if ((this->annealing.bridges > 0) && (this->annealing.state_index != NULL))
  {
    const double all = annealed.evaluate_likelihoods(this->annealing.state_index.get());
    if ((all == -arma::datum::inf) || std::isnan(all))
      return -arma::datum::inf;
  }
  const double now = (this->annealing.bridges > 0) ? annealed.get_cached_likelihoods(index) : annealed.evaluate_likelihoods(index);
  if ((now == -arma::datum::inf) || std::isnan(now))
    return -arma::datum::inf;
  if (this->annealing.bridges == 0)
    log_ratio = now - cached;
  proposed = std::move(annealed);
  return log_ratio;
}

}
