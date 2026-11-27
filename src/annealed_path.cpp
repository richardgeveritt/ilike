#include <algorithm>
#include <cmath>
#include "annealed_path.h"

namespace ilike
{

double log_mean_exp(const std::vector<double> &log_weights)
{
  double top = -arma::datum::inf;
  for (double t : log_weights)
    top = std::max(top, t);
  if (!std::isfinite(top))
    return top;
  double sum = 0.0;
  for (double t : log_weights)
    sum += std::exp(t - top);
  return top + std::log(sum / double(log_weights.size()));
}

FieldPath::FieldPath()
{
}

FieldPath::FieldPath(const std::string &variable_in,
                     UnnormalisedPotentialPtr potential_in,
                     UnnormalisedSweepFunctionPtr sweep_in,
                     UnnormalisedSweepFunctionPtr reverse_sweep_in,
                     const std::vector<std::string> &normaliser_in,
                     size_t bridges_in,
                     size_t sweeps_in,
                     bool reversible_in,
                     Data* data_in)
{
  this->variable = variable_in;
  this->potential = potential_in;
  this->sweep = sweep_in;
  this->reverse_sweep = reverse_sweep_in;
  this->normaliser = normaliser_in;
  this->bridges = bridges_in;
  this->sweeps = std::max(sweeps_in, size_t(1));
  this->reversible = reversible_in;
  this->data = data_in;
  if (this->reversible && (this->reverse_sweep == NULL))
    Rcpp::stop("FieldPath: a reversible path needs the reverse sweep.");
}

AnnealedPath* FieldPath::duplicate() const
{
  return new FieldPath(*this);
}

void FieldPath::set_point(Parameters &state, size_t k, const Parameters &start, const Parameters &end) const
{
  const size_t steps = this->bridges + 1;
  for (const std::string &name : this->normaliser)
  {
    // (the end points exactly: a + (b - a) need not be b in floating point)
    double a = start[name][0], b = end[name][0];
    double value = (k == 0) ? a : ((k == steps) ? b : a + (double(k) / double(steps)) * (b - a));
    state[name] = arma::mat(1, 1, arma::fill::value(value));
  }
}

void FieldPath::move(RandomNumberGenerator &rng, Parameters &state) const
{
  for (size_t s = 0; s < this->sweeps; ++s)
  {
    this->sweep(rng, state, *this->data);
    if (this->reversible)
      this->reverse_sweep(rng, state, *this->data);
  }
}

double FieldPath::single(RandomNumberGenerator &rng, Parameters &state,
                         const Parameters &start, const Parameters &end) const
{
  double log_weight = 0.0;
  for (size_t k = 0; k <= this->bridges; ++k)
  {
    if (k > 0)
    {
      // u_k: the kernel at vartheta_k, from u_{k-1}
      this->set_point(state, k, start, end);
      this->move(rng, state);
    }
    this->set_point(state, k, start, end);
    double here = this->potential(state[this->variable], state, *this->data);
    this->set_point(state, k + 1, start, end);
    log_weight += this->potential(state[this->variable], state, *this->data) - here;
  }
  return log_weight;
}

double FieldPath::population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                             const Parameters &start, const Parameters &end,
                             double threshold) const
{
  const size_t P = states.size();
  if (threshold < 1.0)
  {
    // resampling when the effective sample size falls below threshold x P:
    // the weights are carried between steps, and each step adds the log of
    // the weighted mean incremental weight
    const size_t steps = this->bridges + 1;
    double log_z = 0.0;
    arma::colvec log_w(P, arma::fill::value(-std::log(double(P)))); // (normalised)
    for (size_t t = 1; t <= steps; ++t)
    {
      if (t > 1)
        for (size_t i = 0; i < P; ++i)
        {
          this->set_point(states[i], t - 1, start, end);
          this->move(rng, states[i]);
        }
      arma::colvec incremental(P);
      for (size_t i = 0; i < P; ++i)
      {
        this->set_point(states[i], t - 1, start, end);
        double here = this->potential(states[i][this->variable], states[i], *this->data);
        this->set_point(states[i], t, start, end);
        incremental[i] = this->potential(states[i][this->variable], states[i], *this->data) - here;
      }
      arma::colvec updated = log_w + incremental;
      double top = updated.max();
      if (!std::isfinite(top))
        return top;
      double log_sum = top + std::log(arma::accu(arma::exp(updated - top)));
      log_z += log_sum;
      log_w = updated - log_sum;
      if (t == steps)
        break;
      double ess = 1.0 / arma::accu(arma::exp(2.0 * log_w));
      if (ess < threshold * double(P))
      {
        arma::colvec cumulative = arma::cumsum(arma::exp(log_w));
        std::vector<Parameters> resampled(P);
        for (size_t i = 0; i < P; ++i)
        {
          double u = runif(rng);
          size_t k = std::min(size_t(std::lower_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()), P - 1);
          resampled[i] = states[k].deep_copy();
        }
        states.swap(resampled);
        log_w.fill(-std::log(double(P)));
      }
    }
    return log_z;
  }
  const size_t steps = this->bridges + 1;
  double log_z = 0.0;
  std::vector<double> log_weights(P);
  for (size_t t = 1; t <= steps; ++t)
  {
    for (size_t i = 0; i < P; ++i)
    {
      this->set_point(states[i], t - 1, start, end);
      double here = this->potential(states[i][this->variable], states[i], *this->data);
      this->set_point(states[i], t, start, end);
      log_weights[i] = this->potential(states[i][this->variable], states[i], *this->data) - here;
    }
    log_z += log_mean_exp(log_weights);
    if (t == steps)
      break; // (no move is needed after the last reweighting)
    // multinomial resampling, then the kernel at the new point
    double top = *std::max_element(log_weights.begin(), log_weights.end());
    arma::colvec w(P);
    for (size_t i = 0; i < P; ++i)
      w[i] = std::exp(log_weights[i] - top);
    arma::colvec cumulative = arma::cumsum(w / arma::accu(w));
    std::vector<Parameters> resampled(P);
    for (size_t i = 0; i < P; ++i)
    {
      double u = runif(rng);
      size_t k = std::min(size_t(std::lower_bound(cumulative.begin(), cumulative.end(), u) - cumulative.begin()), P - 1);
      resampled[i] = states[k].deep_copy();
      this->set_point(resampled[i], t, start, end);
      this->move(rng, resampled[i]);
    }
    states.swap(resampled);
  }
  return log_z;
}

}
