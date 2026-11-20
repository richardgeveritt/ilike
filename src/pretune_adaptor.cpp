#include "pretune_adaptor.h"
#include "hmc_proposal_kernel.h"
#include "particle.h"

namespace ilike
{

PretuneAdaptor::PretuneAdaptor()
:ESJDAdaptor(), max_step(0.1), max_steps(100), target(0.9),
 raise_share(0.5), raise_above(0.8), lower_share(0.9), lower_below(0.5), steps_change(5)
{
}

PretuneAdaptor::PretuneAdaptor(const std::string &variable_in, double first_max_step_in, size_t first_max_steps_in,
                               double target_in, bool full_metric_in,
                               double raise_share_in, double raise_above_in,
                               double lower_share_in, double lower_below_in, size_t steps_change_in)
:ESJDAdaptor(variable_in, {"step", "steps"}, {0.0, 1.0}, {first_max_step_in, double(first_max_steps_in)}, {false, true},
             0.0, 0.0, full_metric_in, 1),
 max_step(first_max_step_in), max_steps(first_max_steps_in), target(target_in),
 raise_share(raise_share_in), raise_above(raise_above_in), lower_share(lower_share_in), lower_below(lower_below_in),
 steps_change(steps_change_in)
{
}

PretuneAdaptor::~PretuneAdaptor()
{
}

PretuneAdaptor::PretuneAdaptor(const PretuneAdaptor &another)
:ESJDAdaptor(another)
{
  this->make_copy(another);
}

void PretuneAdaptor::operator=(const PretuneAdaptor &another)
{
  if (this == &another)
    return;
  ESJDAdaptor::operator=(another);
  this->make_copy(another);
}

SMCAdaptor* PretuneAdaptor::duplicate() const
{
  return new PretuneAdaptor(*this);
}

void PretuneAdaptor::make_copy(const PretuneAdaptor &another)
{
  this->max_step = another.max_step;
  this->max_steps = another.max_steps;
  this->target = another.target;
  this->raise_share = another.raise_share;
  this->raise_above = another.raise_above;
  this->lower_share = another.lower_share;
  this->lower_below = another.lower_below;
  this->steps_change = another.steps_change;
}

void PretuneAdaptor::update(RandomNumberGenerator &rng, SMCOutput* current_state, size_t number_of_particles)
{
  // (the values come from the trial moves)
  this->updates = this->updates + 1;
  this->values.set_size(2, number_of_particles);
  this->scores.set_size(number_of_particles);
  this->scores.fill(arma::datum::nan);
  this->compute_metric(current_state);
}

void PretuneAdaptor::median_regression(const arma::colvec &x, const arma::colvec &y, double &a0, double &a1)
{
  arma::uword n = x.n_elem;
  arma::mat X(n, 2);
  X.col(0).ones();
  X.col(1) = x;
  arma::colvec w = arma::ones<arma::colvec>(n);
  arma::colvec beta = arma::solve(X, y);
  double scale = std::max(1e-12, arma::median(arma::abs(y)));
  for (int k = 0; k < 100; ++k)
  {
    arma::colvec r = arma::abs(y - X*beta);
    w = 1.0/arma::max(r, arma::colvec(n, arma::fill::value(1e-8*scale)));
    arma::mat Xw = X.each_col() % arma::sqrt(w);
    arma::colvec next;
    if (!arma::solve(next, Xw, y % arma::sqrt(w)))
      break;
    bool converged = arma::norm(next - beta) <= 1e-10*(1.0 + arma::norm(beta));
    beta = next;
    if (converged)
      break;
  }
  a0 = beta[0];
  a1 = beta[1];
}

void PretuneAdaptor::trial_moves(RandomNumberGenerator &rng, size_t N, const std::function<Particle(size_t)> &start)
{
  const HMCProposalKernel* kernel = dynamic_cast<const HMCProposalKernel*>(this->proposal);
  if (kernel == NULL)
    Rcpp::stop("pretune tunes ilike::hmc (in metropolis_hastings).");
  arma::colvec eps(N), steps(N), energy(N), scores_trial(N);
  double trial_values[2];
  ParticleTuning trial;
  ParticleTuning::Row row;
  row.instance = kernel->get_instance_index();
  row.values = trial_values;
  row.adaptor = this;
  trial.rows.push_back(row);
  for (size_t i = 0; i < N; ++i)
  {
    double e = runif(rng)*this->max_step;
    if (!(e > 0.0))
      e = 1e-3*this->max_step;
    double L = std::floor(runif(rng)*double(this->max_steps)) + 1.0;
    if (L > double(this->max_steps))
      L = double(this->max_steps);
    eps[i] = e;
    steps[i] = L;
    trial_values[0] = e;
    trial_values[1] = L;
    Particle particle = start(i);
    particle.tuning = &trial;
    const Parameters &base = particle.get_transformed_parameters(kernel);
    arma::colvec x0 = base.get_colvec(kernel->block());
    double log_kinetic_ratio;
    Parameters proposed = kernel->trajectory(rng, particle, log_kinetic_ratio);
    arma::colvec x1 = proposed.get_colvec(kernel->block());
    double dE = -arma::datum::inf;
    if (std::isfinite(log_kinetic_ratio))
    {
      double L0 = particle.get_cached_likelihoods(kernel->get_index()) + kernel->log_jacobian(base);
      double L1 = kernel->target_at(particle, base, x1);
      dE = L1 - L0 + log_kinetic_ratio;
      if (std::isnan(dE))
        dE = -arma::datum::inf;
    }
    energy[i] = dE;
    double alpha = std::isfinite(dE) ? std::min(1.0, std::exp(dE)) : 0.0;
    double jump = (alpha > 0.0) ? this->mahalanobis(x1 - x0) : 0.0;
    scores_trial[i] = std::isfinite(jump) ? alpha*jump/L : 0.0;
    particle.tuning = NULL;
  }

  // eps*: the median regression of |dE| on eps^2 over the finite energy errors
  arma::uvec finite = arma::find_finite(energy);
  double wanted = std::abs(std::log(this->target));
  double next_max = this->max_step;
  if (finite.n_elem >= 10)
  {
    double a0, a1;
    median_regression(arma::square(eps.elem(finite)), arma::abs(energy.elem(finite)), a0, a1);
    if ((a1 > 0.0) && (wanted > a0))
      next_max = std::sqrt((wanted - a0)/a1);
    else if (a0 >= wanted)
      next_max = 0.5*this->max_step;
    else
      next_max = 2.0*this->max_step;
  }
  else
    next_max = 0.5*this->max_step; // (most trial trajectories failed)

  // the real move's values: the trial pairs resampled with weights Lambda
  double total = arma::accu(scores_trial);
  arma::colvec weights = (total > 0.0) && std::isfinite(total) ? arma::colvec(scores_trial/total)
                                                              : arma::colvec(N, arma::fill::value(1.0/double(N)));
  arma::colvec cumulative = arma::cumsum(weights);
  std::vector<size_t> chosen(N);
  double u = runif(rng);
  size_t j = 0;
  for (size_t i = 0; i < N; ++i)
  {
    double point = (double(i) + u)/double(N);
    while ((j + 1 < N) && (cumulative[j] < point))
      ++j;
    chosen[i] = j;
  }
  for (size_t i = N; i > 1; --i)
  {
    size_t k = std::min(size_t(std::floor(runif(rng)*double(i))), i - 1);
    std::swap(chosen[i - 1], chosen[k]);
  }
  size_t high = 0, low = 0;
  for (size_t i = 0; i < N; ++i)
  {
    this->values(0, i) = eps[chosen[i]];
    this->values(1, i) = steps[chosen[i]];
    if (steps[chosen[i]] > this->raise_above*double(this->max_steps))
      ++high;
    if (steps[chosen[i]] < this->lower_below*double(this->max_steps))
      ++low;
  }
  if (double(high) > this->raise_share*double(N))
    this->max_steps = this->max_steps + this->steps_change;
  else if (double(low) > this->lower_share*double(N))
    this->max_steps = std::max<size_t>(std::max<size_t>(1, this->steps_change),
                                       this->max_steps - std::min(this->max_steps, this->steps_change));
  this->max_step = next_max;
}

void PretuneAdaptor::summary(std::vector<std::string> &names_out, std::vector<double> &numbers) const
{
  ESJDAdaptor::summary(names_out, numbers);
  // (the ranges for the next step's trial moves)
  names_out.push_back("step_max");
  numbers.push_back(this->max_step);
  names_out.push_back("steps_max");
  numbers.push_back(double(this->max_steps));
}

}
