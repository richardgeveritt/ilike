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
  this->resampling_threshold = another.resampling_threshold;
  this->log_normaliser_at = another.log_normaliser_at;
  this->replicates = another.replicates;
  this->tempering_variable = another.tempering_variable;
  this->engine_single = another.engine_single;
  this->engine_population = another.engine_population;
  this->engine_normalise = another.engine_normalise;
  this->plugin = another.plugin;
  this->plugin_particles = another.plugin_particles;
  this->plugin_steps = another.plugin_steps;
  this->plugin_sweeps = another.plugin_sweeps;
  this->plugin_resampling = another.plugin_resampling;
  this->roulette = another.roulette;
  this->roulette_continue = another.roulette_continue;
  this->roulette_scale = another.roulette_scale;
  this->block_poisson = another.block_poisson;
  this->block_poisson_blocks = another.block_poisson_blocks;
  this->block_poisson_mean = another.block_poisson_mean;
  this->block_poisson_lower = another.block_poisson_lower;
  this->block_poisson_correlated = another.block_poisson_correlated;
}

double UnnormalisedNormaliserEstimator::draw_replicate_log_estimate(const Parameters &parameters, double *sign) const
{
  if (sign != NULL)
    *sign = 1.0;
  // (dev/plan.md 3.22 U9 R1) the plug-in: -log Zhat(theta), from the base measure
  if (this->plugin)
  {
    Parameters end;
    for (const std::string &name : this->normaliser)
      end[name] = parameters[name];
    if (!this->roulette)
      return -this->base_log_evidence(*this->rng, end, this->plugin_particles, this->plugin_steps, this->plugin_sweeps,
                                      this->plugin_resampling);
    // (dev/plan.md 3.22 U9 R2) Russian roulette: the scale from a pilot
    // estimate, then the randomly truncated series
    auto log_z = [&]() {
      return this->base_log_evidence(*this->rng, end, this->plugin_particles, this->plugin_steps, this->plugin_sweeps,
                                     this->plugin_resampling);
    };
    const double log_c = -std::log(this->roulette_scale) - log_z();
    double total = 1.0, product = 1.0, survival = 1.0;
    for (size_t n = 1; n <= 100000; ++n)
    {
      if (!(runif(*this->rng) < this->roulette_continue))
        break;
      product = product * (1.0 - std::exp(log_c + log_z()));
      survival = survival * this->roulette_continue;
      total = total + product / survival;
    }
    if (sign != NULL)
      *sign = (total < 0.0) ? -1.0 : 1.0;
    return log_c + std::log(std::abs(total));
  }
  // (dev/plan.md 3.22 U3c) each point an annealed importance sampling run
  // along the path to theta_hat (with no bridges, gamma(u | theta_hat)/gamma(u | theta));
  // with resampling, the runs moved together (an SMC)
  if (this->resampling_threshold > 0.0)
  {
    // the M runs moved together, with resampling (dev/plan.md 3.22 U3c)
    std::vector<Parameters> states(this->points);
    for (size_t m = 0; m < this->points; ++m)
    {
      arma::mat field = this->draw(parameters);
      if (field.has_nan())
        return -arma::datum::inf;
      states[m] = parameters.deep_copy();
      states[m][this->variable] = field;
    }
    // (resampling when the effective sample size falls below half)
    if (this->engine_population)
      return this->engine_population->population(*this->rng, states, parameters, this->at_parameters()) - this->log_normaliser_at;
    return this->path().population(*this->rng, states, parameters, this->at_parameters(), this->resampling_threshold) - this->log_normaliser_at;
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

void UnnormalisedNormaliserEstimator::set_normalise(size_t particles_in, size_t steps_in, size_t sweeps_in, size_t seed_value,
                                                    double resampling_threshold_in)
{
  RandomNumberGenerator rng;
  rng.seed(seed_value ^ 0x9e3779b97f4a7c15ULL);
  this->log_normaliser_at = this->base_log_evidence(rng, this->at_parameters(), particles_in, steps_in, sweeps_in, resampling_threshold_in);
}

void UnnormalisedNormaliserEstimator::set_roulette(double continue_in, double scale_in)
{
  this->roulette = true;
  this->roulette_continue = continue_in;
  this->roulette_scale = scale_in;
}

void UnnormalisedNormaliserEstimator::set_block_poisson(size_t blocks_in, double mean_in, double lower_in, bool correlated_in)
{
  if (this->replicates > 1)
    Rcpp::stop("debias = block_poisson(…) with replicates is not available yet.");
  this->block_poisson = true;
  this->block_poisson_blocks = std::max(blocks_in, size_t(1));
  this->block_poisson_mean = mean_in;
  this->block_poisson_lower = lower_in;
  this->block_poisson_correlated = correlated_in;
}

double UnnormalisedNormaliserEstimator::plugin_log_z(RandomNumberGenerator &rng, const Parameters &parameters) const
{
  Parameters end;
  for (const std::string &name : this->normaliser)
    end[name] = parameters[name];
  return this->base_log_evidence(rng, end, this->plugin_particles, this->plugin_steps, this->plugin_sweeps,
                                 this->plugin_resampling);
}

void UnnormalisedNormaliserEstimator::set_plugin(size_t particles_in, size_t steps_in, size_t sweeps_in, double resampling_threshold_in)
{
  this->plugin = true;
  this->plugin_particles = particles_in;
  this->plugin_steps = steps_in;
  this->plugin_sweeps = sweeps_in;
  this->plugin_resampling = resampling_threshold_in;
}

double UnnormalisedNormaliserEstimator::base_log_evidence(RandomNumberGenerator &rng, const Parameters &end, size_t particles_in,
                                                          size_t steps_in, size_t sweeps_in, double resampling_threshold_in) const
{
  if (this->engine_normalise)
  {
    // the general engine: an SMC from the base measure to `end`
    return this->engine_normalise->log_evidence(rng, end);
  }
  const size_t P = std::max(particles_in, size_t(2));
  const size_t T = std::max(steps_in, size_t(1));
  // the path from 0 to `end` in the normaliser's parameters
  Parameters start;
  for (size_t j = 0; j < this->normaliser.size(); ++j)
    start[this->normaliser[j]] = arma::mat(arma::size(end[this->normaliser[j]]), arma::fill::zeros);
  // P draws from the base measure: a sweep at beta = 0 (from an empty field,
  // the generated sweep first fills it from the base measure)
  std::vector<Parameters> states(P);
  for (size_t i = 0; i < P; ++i)
  {
    states[i] = start.deep_copy();
    this->sweep(rng, states[i], *this->current_data);
  }
  // an SMC with resampling along the path (the shared engine, dev/plan.md 3.22 U3c)
  FieldPath path(this->variable, this->potential, this->sweep, NULL, this->normaliser, T - 1, sweeps_in, false, this->current_data);
  return path.population(rng, states, start, end, resampling_threshold_in);
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
  if (this->engine_single || this->engine_population)
  {
    // the general engine: its kernel at the parameters, from the latent field
    // or from a draw from the base measure
    const SMCPath *engine = this->engine_single ? this->engine_single.get() : this->engine_population.get();
    Parameters state = parameters.deep_copy();
    if (state.find_vector(this->variable) == NULL)
      state[this->variable] = engine->draw_base(*this->rng, parameters)[this->variable];
    engine->draw(*this->rng, state, parameters, std::max(this->sweeps, size_t(1)));
    return state[this->variable];
  }
  // sweeps (approximate): from the field in the parameters if it is there (a
  // latent field), otherwise from the base measure
  Parameters state = parameters.deep_copy();
  for (size_t k = 0; k < std::max(this->sweeps, size_t(1)); ++k)
    this->sweep(*this->rng, state, *this->current_data);
  return state[this->variable];
}


double UnnormalisedNormaliserEstimator::log_mean_exp(const std::vector<double> &log_weights)
{
  return ilike::log_mean_exp(log_weights);
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
  if (this->engine_single)
    return this->engine_single->single(*this->rng, state, parameters, this->at_parameters());
  return this->path().single(*this->rng, state, parameters, this->at_parameters());
}

FieldPath UnnormalisedNormaliserEstimator::path() const
{
  return FieldPath(this->variable, this->potential, this->sweep, NULL, this->normaliser, this->bridges, this->bridge_sweeps, false, this->current_data);
}

Parameters UnnormalisedNormaliserEstimator::at_parameters() const
{
  Parameters end;
  for (size_t j = 0; j < this->normaliser.size(); ++j)
    end[this->normaliser[j]] = arma::mat(1, 1, arma::fill::value(this->at[j]));
  return end;
}

void UnnormalisedNormaliserEstimator::set_engine(SMCPath* single_in, SMCPath* population_in, SMCPath* normalise_in, size_t draw_sweeps_in)
{
  this->engine_single.reset(single_in);
  this->engine_population.reset(population_in);
  this->engine_normalise.reset(normalise_in);
  this->sweeps = std::max(draw_sweeps_in, size_t(1));
}

}
