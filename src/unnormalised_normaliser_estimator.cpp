#include <algorithm>
#include "unnormalised_normaliser_estimator.h"
#include "unnormalised_normaliser_estimator_output.h"

namespace ilike
{

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator()
  :LikelihoodEstimator()
{
  this->potential = NULL;
  this->exact = NULL;
  this->sweep = NULL;
  this->sweeps = 0;
  this->points = 1;
}

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator(RandomNumberGenerator* rng_in,
                                                                 size_t* seed_in,
                                                                 Data* data_in,
                                                                 const std::string &variable_in,
                                                                 UnnormalisedPotentialPtr potential_in,
                                                                 UnnormalisedExactFunctionPtr exact_in,
                                                                 UnnormalisedSweepFunctionPtr sweep_in,
                                                                 size_t sweeps_in,
                                                                 size_t points_in,
                                                                 const std::vector<std::string> &normaliser_in,
                                                                 const std::vector<double> &at_in,
                                                                 bool smcfixed_flag_in)
  :LikelihoodEstimator(rng_in, seed_in, data_in, Parameters(), smcfixed_flag_in)
{
  this->variable = variable_in;
  this->potential = potential_in;
  this->exact = exact_in;
  this->sweep = sweep_in;
  this->sweeps = sweeps_in;
  this->points = std::max(points_in, size_t(1));
  this->normaliser = normaliser_in;
  this->at = at_in;
}

UnnormalisedNormaliserEstimator::~UnnormalisedNormaliserEstimator()
{
}

UnnormalisedNormaliserEstimator::UnnormalisedNormaliserEstimator(const UnnormalisedNormaliserEstimator &another)
  :LikelihoodEstimator(another)
{
  this->make_copy(another);
}

void UnnormalisedNormaliserEstimator::operator=(const UnnormalisedNormaliserEstimator &another)
{
  if (this == &another)
    return;
  LikelihoodEstimator::operator=(another);
  this->make_copy(another);
}

LikelihoodEstimator* UnnormalisedNormaliserEstimator::duplicate() const
{
  return new UnnormalisedNormaliserEstimator(*this);
}

void UnnormalisedNormaliserEstimator::make_copy(const UnnormalisedNormaliserEstimator &another)
{
  this->variable = another.variable;
  this->potential = another.potential;
  this->exact = another.exact;
  this->sweep = another.sweep;
  this->sweeps = another.sweeps;
  this->points = another.points;
  this->normaliser = another.normaliser;
  this->at = another.at;
  this->bridges = another.bridges;
  this->bridge_sweeps = another.bridge_sweeps;
  this->log_normaliser_at = another.log_normaliser_at;
  this->replicates = another.replicates;
  this->tempering_variable = another.tempering_variable;
}

double UnnormalisedNormaliserEstimator::draw_replicate_log_estimate(const Parameters &parameters) const
{
  if (this->bridges == 0)
  {
    std::vector<arma::mat> fields(this->points);
    for (size_t m = 0; m < this->points; ++m)
      fields[m] = this->draw(parameters);
    return this->log_estimate(fields, parameters);
  }
  std::vector<double> log_weights(this->points);
  for (size_t m = 0; m < this->points; ++m)
    log_weights[m] = this->annealed_log_weight(this->draw(parameters), parameters);
  return log_mean_exp(log_weights) - this->log_normaliser_at;
}

size_t UnnormalisedNormaliserEstimator::active_replicates(const Parameters &parameters) const
{
  if (this->tempering_variable.empty())
    return this->replicates;
  const std::vector<std::string> names = parameters.get_vector_variables();
  if (std::find(names.begin(), names.end(), this->tempering_variable) == names.end())
    return this->replicates;
  double whole = parameters[this->tempering_variable][0] * double(this->replicates);
  double k = std::round(whole);
  if (std::abs(whole - k) > 1e-8)
    Rcpp::stop("normaliser = sav(…) with replicates under data tempering: the sequence must add whole replicates (data(…, anneal = false) "
               "with a schedule of multiples of 1/replicates); the estimate of a fraction of a replicate's normalising constant is not available.");
  return std::min(size_t(std::max(k, 0.0)), this->replicates);
}

void UnnormalisedNormaliserEstimator::set_normalise(size_t particles_in, size_t steps_in, size_t sweeps_in, size_t seed_value)
{
  const size_t P = std::max(particles_in, size_t(2));
  const size_t T = std::max(steps_in, size_t(1));
  const size_t m = std::max(sweeps_in, size_t(1));
  RandomNumberGenerator rng;
  rng.seed(seed_value ^ 0x9e3779b97f4a7c15ULL);
  // the parameters at a point of the path (only the normaliser's parameters
  // are read by the field's functions)
  auto at_point = [&](double s) {
    Parameters point;
    for (size_t j = 0; j < this->normaliser.size(); ++j)
      point[this->normaliser[j]] = arma::mat(1, 1, arma::fill::value(s * this->at[j]));
    return point;
  };
  // P draws from the base measure: a sweep at beta = 0 (from an empty field,
  // the generated sweep first fills it from the base measure)
  std::vector<Parameters> particles(P);
  Parameters start = at_point(0.0);
  for (size_t i = 0; i < P; ++i)
  {
    particles[i] = start.deep_copy();
    this->sweep(rng, particles[i], *this->current_data);
  }
  double log_z = 0.0;
  std::vector<double> log_weights(P);
  for (size_t t = 1; t <= T; ++t)
  {
    Parameters previous = at_point(double(t - 1) / double(T)), current = at_point(double(t) / double(T));
    for (size_t i = 0; i < P; ++i)
    {
      const arma::mat &field = particles[i][this->variable];
      log_weights[i] = this->potential(field, current, *this->current_data) - this->potential(field, previous, *this->current_data);
    }
    log_z += log_mean_exp(log_weights);
    // multinomial resampling, then m sweeps at the new point
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
      resampled[i] = particles[k].deep_copy();
      resampled[i].deep_overwrite_with_variables_in_argument(current);
      for (size_t s = 0; s < m; ++s)
        this->sweep(rng, resampled[i], *this->current_data);
    }
    particles.swap(resampled);
  }
  this->log_normaliser_at = log_z;
}

void UnnormalisedNormaliserEstimator::set_bridges(size_t bridges_in, size_t bridge_sweeps_in)
{
  this->bridges = bridges_in;
  this->bridge_sweeps = std::max(bridge_sweeps_in, size_t(1));
}

LikelihoodEstimatorOutput* UnnormalisedNormaliserEstimator::initialise()
{
  return new UnnormalisedNormaliserEstimatorOutput(this);
}

LikelihoodEstimatorOutput* UnnormalisedNormaliserEstimator::initialise(const Parameters &parameters)
{
  return this->initialise();
}

void UnnormalisedNormaliserEstimator::setup()
{
}

void UnnormalisedNormaliserEstimator::setup(const Parameters &parameters)
{
}

void UnnormalisedNormaliserEstimator::specific_change_data(Data* new_data)
{
}

arma::mat UnnormalisedNormaliserEstimator::draw(const Parameters &parameters) const
{
  if (this->exact != NULL)
  {
    Parameters drawn = this->exact(*this->rng, parameters, *this->current_data);
    return drawn[this->variable];
  }
  // sweeps (approximate): from the field in the parameters if it is there (a
  // latent field), otherwise from the base measure
  Parameters state = parameters.deep_copy();
  for (size_t k = 0; k < std::max(this->sweeps, size_t(1)); ++k)
    this->sweep(*this->rng, state, *this->current_data);
  return state[this->variable];
}

double UnnormalisedNormaliserEstimator::log_estimate(const std::vector<arma::mat> &fields, const Parameters &parameters) const
{
  Parameters at_parameters = parameters.deep_copy();
  for (size_t k = 0; k < this->normaliser.size(); ++k)
    at_parameters[this->normaliser[k]] = arma::mat(1, 1, arma::fill::value(this->at[k]));
  // log-sum-exp of log gamma(u | theta_hat) - log gamma(u | theta)
  // (a draw that could not be made, e.g. exact draws of an ising field at
  // beta < 0: density 0 at these parameters)
  for (const arma::mat &field : fields)
    if (field.has_nan())
      return -arma::datum::inf;
  std::vector<double> terms(fields.size());
  for (size_t m = 0; m < fields.size(); ++m)
    terms[m] = this->potential(fields[m], at_parameters, *this->current_data) -
               this->potential(fields[m], parameters, *this->current_data);
  return log_mean_exp(terms) - this->log_normaliser_at;
}

double UnnormalisedNormaliserEstimator::log_mean_exp(const std::vector<double> &log_weights)
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

// Annealed importance sampling (Neal 2001) from the field at theta to the
// auxiliary distribution at theta_hat: with vartheta_k = theta + k/(K+1)
// (theta_hat - theta), k = 0, ..., K+1, u_0 = field (drawn at theta) and u_k
// made from u_{k-1} by sweeps that leave the field at vartheta_k invariant,
// the weight sum_{k=0}^{K} log gamma(u_k | vartheta_{k+1}) - log gamma(u_k | vartheta_k)
// has exp(weight) unbiased for Z(theta_hat)/Z(theta). (For ising and potts,
// whose potentials are linear in beta, this straight path is the geometric
// bridge gamma(.|theta)^(1-s) gamma(.|theta_hat)^s.)
double UnnormalisedNormaliserEstimator::annealed_log_weight(const arma::mat &field, const Parameters &parameters) const
{
  if (field.has_nan())
    return -arma::datum::inf;
  Parameters state = parameters.deep_copy();
  state[this->variable] = field;
  const double steps = double(this->bridges + 1);
  auto set_point = [&](size_t k) {
    for (size_t j = 0; j < this->normaliser.size(); ++j)
    {
      double start = parameters[this->normaliser[j]][0];
      state[this->normaliser[j]] = arma::mat(1, 1, arma::fill::value(start + (double(k) / steps) * (this->at[j] - start)));
    }
  };
  double log_weight = 0.0;
  for (size_t k = 0; k <= this->bridges; ++k)
  {
    if (k > 0)
    {
      // u_k: sweeps at vartheta_k, from u_{k-1}
      set_point(k);
      for (size_t s = 0; s < this->bridge_sweeps; ++s)
        this->sweep(*this->rng, state, *this->current_data);
    }
    set_point(k);
    double here = this->potential(state[this->variable], state, *this->current_data);
    set_point(k + 1);
    log_weight += this->potential(state[this->variable], state, *this->current_data) - here;
  }
  return log_weight;
}

}
