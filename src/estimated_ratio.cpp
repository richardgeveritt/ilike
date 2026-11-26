#include <cmath>
#include "estimated_ratio.h"
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
  if ((this->step != NULL) && this->step->changes_state())
    Rcpp::stop("EstimatedRatio: auxiliary draws that change the state (a latent z', reversible jump) not written yet.");
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
