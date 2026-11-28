#include <cmath>
#include "estimated_ratio.h"
#include "annealed_path.h"
#include "particle.h"
#include "chunk_parallel.h"
#include "unnormalised_normaliser_estimator.h"

namespace ilike
{

EstimatedRatio::EstimatedRatio()
{
}

EstimatedRatio::EstimatedRatio(RatioStep* step_in)
{
  this->step.reset(step_in);
}

// (dev/plan.md 3.22 U6, C313) the averaged form for state steps with a
// selection among the N new states (RmJ-MCMC, 2101.01253 Algorithm 3, with
// the direction probabilities omega of its Remark):
//  c = 1 (probability omega(theta -> vartheta, 1)): N draws to vartheta, log
//    targets t_i; the new state the i-th with probability proportional to
//    exp(t_i); log term = log omega(vartheta -> theta, 2) - log omega(theta -> vartheta, 1)
//    + log (1/N) sum_i exp(t_i);
//  c = 2: one draw to vartheta, new state z', t_1 = -(its log target) (the
//    exact reverse's), and N - 1 draws from (vartheta, z') back to theta,
//    t_i; log term = -[log omega(theta -> vartheta, 2) - log omega(vartheta -> theta, 1)
//    + log (1/N) sum_i exp(t_i)].
// The proposal's density ratio is the same for every draw (it reads only
// theta and vartheta), so the move adds it outside the average.
static Particle back_towards(const Particle &from_state, const Particle &current, const Particle &proposed)
{
  // from_state with the variables the proposal changed set back to their current values
  Particle shell = from_state;
  const Parameters &before = current.parameters;
  const Parameters &after = proposed.parameters;
  for (auto a = after.vector_begin(); a != after.vector_end(); ++a)
  {
    const arma::mat* b = before.find_vector(a->first);
    if ((b == NULL) || (b == a->second.first.get()) || ((b->n_elem == a->second.first->n_elem) && arma::approx_equal(*b, *a->second.first, "absdiff", 0.0)))
      continue;
    shell.parameters.overwrite_variable(a->first, *b);
  }
  return shell;
}

double EstimatedRatio::averaged_state_log_term(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const
{
  if (!this->step->averages())
    Rcpp::stop("points = N is not available with this move (its draws change the state, and it has no averaged form).");
  const size_t N = this->points;
  const double omega_forward_1 = this->step->first_direction_probability(current, proposed);
  const bool first = runif(rng) < omega_forward_1;
  std::vector<double> terms(N);
  if (first)
  {
    std::vector<Particle> states(N, proposed);
    for (size_t i = 0; i < N; ++i)
      terms[i] = this->step->draw_state(rng, index, current, states[i]);
    // the selection, with probabilities proportional to exp(t_i)
    double largest = -arma::datum::inf;
    for (double t : terms)
      largest = std::max(largest, t);
    if (!std::isfinite(largest))
      return -arma::datum::inf;
    double total = 0.0;
    for (double t : terms)
      total += std::exp(t - largest);
    double u = runif(rng) * total, cumulative = 0.0;
    size_t chosen = N - 1;
    for (size_t i = 0; i < N; ++i)
    {
      cumulative += std::exp(terms[i] - largest);
      if (u < cumulative)
      {
        chosen = i;
        break;
      }
    }
    proposed = std::move(states[chosen]);
    Particle shell = back_towards(proposed, current, proposed);
    const double omega_back_2 = 1.0 - this->step->first_direction_probability(proposed, shell);
    return std::log(omega_back_2) - std::log(omega_forward_1) + ilike::log_mean_exp(terms);
  }
  Particle forward = proposed;
  const double t_forward = this->step->draw_state(rng, index, current, forward);
  if (!std::isfinite(t_forward))
    return -arma::datum::inf;
  terms[0] = -t_forward;
  Particle shell = back_towards(forward, current, proposed);
  for (size_t i = 1; i < N; ++i)
  {
    Particle back = shell;
    terms[i] = this->step->draw_state(rng, index, forward, back);
  }
  const double omega_back_1 = this->step->first_direction_probability(forward, shell);
  proposed = std::move(forward);
  return -(std::log(1.0 - omega_forward_1) - std::log(omega_back_1) + ilike::log_mean_exp(terms));
}

double EstimatedRatio::state_log_term(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed,
                                      Particle* refreshed, bool* refreshed_made) const
{
  if (refreshed_made != NULL)
    *refreshed_made = false;
  if (this->points != 1)
    return this->averaged_state_log_term(rng, index, current, proposed);
  if (!this->step->uses_direction() && !this->step->refreshes())
    return this->step->draw_state(rng, index, current, proposed);
  // (C312) the direction, for a step that uses it
  int direction = 1;
  if (this->step->uses_direction())
    direction = (runif(rng) < 0.5) ? 1 : 2;
  return this->step->draw_state(rng, index, current, proposed, direction, refreshed, refreshed_made);
}

EstimatedRatio::EstimatedRatio(const EstimatedRatio &another)
{
  *this = another;
}

EstimatedRatio& EstimatedRatio::operator=(const EstimatedRatio &another)
{
  if (this == &another)
    return *this;
  if (another.step != NULL)
    this->step.reset(another.step->duplicate());
  else
    this->step.reset();
  this->points = another.points;
  this->noisy = another.noisy;
  this->parallel = another.parallel;
  this->number_of_threads = another.number_of_threads;
  return *this;
}

void EstimatedRatio::set_parallel(bool parallel_in, size_t number_of_threads_in)
{
  this->parallel = parallel_in;
  this->number_of_threads = std::max(number_of_threads_in, size_t(1));
}

double EstimatedRatio::log_term(RandomNumberGenerator &rng,
                                const Parameters &current,
                                const Parameters &proposed) const
{
  // one draw: the step's term (for exchange, the exchange algorithm)
  if (this->points == 1)
    return this->step->draw(rng, current, proposed);
  // N draws: the direction (averaged) or always the first (noisy), then the
  // draws in their own streams
  const bool forward = this->noisy || (runif(rng) < 0.5);
  RandomNumberGenerator first;
  first.seed(rng());
  std::vector<RandomNumberGenerator> streams = chunk_streams(first, this->points);
  std::vector<double> terms(this->points);
  const RatioStep* step = this->step.get();
  std::function<void(size_t, size_t, size_t)> body = [&](size_t g, size_t a, size_t b) {
    for (size_t i = a; i < b; ++i)
    {
      if (forward)
        terms[i] = step->draw(streams[i], current, proposed);          // c = 1
      else if (i == 0)
        terms[i] = -step->draw(streams[i], current, proposed);         // c = 2: the forward draw
      else
        terms[i] = step->draw(streams[i], proposed, current);          // c = 2: backward draws
    }
  };
  run_chunk_groups(0, this->points, this->parallel, this->number_of_threads, body);
  if (forward)
    return UnnormalisedNormaliserEstimator::log_mean_exp(terms);
  if (std::isinf(terms[0]))
    return -arma::datum::inf; // (no draw at the proposed parameters)
  return -UnnormalisedNormaliserEstimator::log_mean_exp(terms);
}

}
