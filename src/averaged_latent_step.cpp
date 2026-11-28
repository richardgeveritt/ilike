#include <cmath>
#include "averaged_latent_step.h"
#include "particle.h"
#include "vector_index.h"
#include "vector_factor_variables.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "utils.h"

namespace ilike
{

AveragedLatentStep::AveragedLatentStep(const std::string &variable_in, size_t points_in, bool refresh_in,
                                       SimulateIndependentProposalPtr simulate_in,
                                       const std::vector<size_t> &latent_factors_in, const std::vector<size_t> &own_factors_in)
{
  this->variable = variable_in;
  this->points = std::max(points_in, size_t(2));
  this->refresh = refresh_in;
  this->simulate = simulate_in;
  this->latent_factors = latent_factors_in;
  this->own_factors = own_factors_in;
  this->latent_index = std::make_shared<VectorIndex>(latent_factors_in);
  this->own_index = std::make_shared<VectorIndex>(own_factors_in);
}

AveragedLatentStep::AveragedLatentStep(const AveragedLatentStep &another)
{
  this->variable = another.variable;
  this->points = another.points;
  this->refresh = another.refresh;
  this->simulate = another.simulate;
  this->latent_factors = another.latent_factors;
  this->own_factors = another.own_factors;
  this->latent_index = another.latent_index;
  this->own_index = another.own_index;
}

AveragedLatentStep::~AveragedLatentStep()
{
}

RatioStep* AveragedLatentStep::duplicate() const
{
  return new AveragedLatentStep(*this);
}

double AveragedLatentStep::draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const
{
  Rcpp::stop("AveragedLatentStep::draw: an averaged move's step changes the state (draw_state).");
}

double AveragedLatentStep::draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const
{
  return this->draw_state(rng, index, current, proposed, 1, NULL, NULL);
}

// (the particle's proposal stores belong to its old values: made again, as a new particle's are)
static void remake_stores(Particle &particle)
{
  particle.current_proposal_store.clear();
  particle.simulate_proposal_variables(particle.proposals_to_transform_for_pointer, particle.proposals_to_find_gradient_for_pointer);
}

void AveragedLatentStep::element_changes(Particle &particle, const std::vector<arma::mat> &candidates, arma::mat &full, arma::mat &own) const
{
  VectorFactorVariables* factor_variables = dynamic_cast<VectorFactorVariables*>(particle.factor_variables.get());
  if (factor_variables == NULL)
    Rcpp::stop("averaged(…, latent = " + this->variable + "): the model's factors are needed as a list of outputs.");
  // every factor reading z made of terms that each read at most one element
  for (size_t f : this->latent_factors)
  {
    LikelihoodEstimatorOutput* output = factor_variables->get_output(f);
    if (output == NULL)
      continue;
    int separate = (output->get_likelihood_estimator() != NULL) ? output->get_likelihood_estimator()->terms_read_one_element(this->variable) : -1;
    if (separate == 0)
      Rcpp::stop("averaged(…, latent = " + this->variable + "): factor " + std::to_string(f + 1) + " has a term reading more than one element of " +
                 this->variable + " (or all of it); the averaged move needs the elements of " + this->variable + " independent given the other variables.");
    if (separate < 0)
      Rcpp::stop("averaged(…, latent = " + this->variable + "): factor " + std::to_string(f + 1) + " reads " + this->variable +
                 " but is not written as model statements, so its terms reading each element are not known.");
  }
  // the latent's matrix owned by this particle alone, so that it can be changed in place
  {
    Parameters mine;
    mine[this->variable] = static_cast<const Parameters&>(particle.parameters)[this->variable];
    particle.parameters.deep_overwrite_with_variables_in_argument(mine);
  }
  arma::mat &x = particle.parameters[this->variable];
  const double base_full = particle.get_cached_likelihoods(this->latent_index.get());
  const double base_own = particle.get_cached_likelihoods(this->own_index.get());
  for (size_t f : this->latent_factors)
  {
    LikelihoodEstimatorOutput* output = factor_variables->get_output(f);
    if (output != NULL)
      output->prepare_in_place(particle.parameters);
  }
  const size_t T = x.n_elem, M = candidates.size();
  full.zeros(T, M);
  own.zeros(T, M);
  std::vector<LikelihoodEstimatorOutput*> changed;
  for (size_t t = 0; t < T; ++t)
  {
    const size_t row = t % x.n_rows, column = t / x.n_rows;
    const double old_value = x(row, column);
    for (size_t i = 1; i < M; ++i)
    {
      x(row, column) = candidates[i](row, column);
      changed.clear();
      for (size_t f : this->latent_factors)
      {
        LikelihoodEstimatorOutput* output = factor_variables->get_output(f);
        if (output == NULL)
          continue;
        if (!output->change_in_place(particle.parameters, this->variable, row, column))
          Rcpp::stop("averaged(…, latent = " + this->variable + "): factor " + std::to_string(f + 1) + " cannot be changed one element at a time.");
        changed.push_back(output);
      }
      full(t, i) = particle.get_cached_likelihoods(this->latent_index.get()) - base_full;
      own(t, i) = particle.get_cached_likelihoods(this->own_index.get()) - base_own;
      for (auto c = changed.rbegin(); c != changed.rend(); ++c)
        (*c)->undo_in_place();
    }
    x(row, column) = old_value;
  }
}

void AveragedLatentStep::install(Particle &particle, const std::vector<arma::mat> &candidates, const std::vector<size_t> &chosen) const
{
  arma::mat z = candidates[0];
  for (size_t t = 0; t < z.n_elem; ++t)
    z(t) = candidates[chosen[t]](t);
  particle.parameters.overwrite_variable(this->variable, z);
  particle.evaluate_likelihoods(this->latent_index.get());
  remake_stores(particle);
}

// a draw from the probabilities proportional to exp(log_weights)
static size_t draw_from(RandomNumberGenerator &rng, const arma::rowvec &log_weights)
{
  const double largest = log_weights.max();
  double u = runif(rng);
  if (!std::isfinite(largest))
    return 0;
  arma::rowvec w = arma::exp(log_weights - largest);
  w = w / arma::accu(w);
  double cumulative = 0.0;
  for (size_t i = 0; i < w.n_elem; ++i)
  {
    cumulative += w(i);
    if (u < cumulative)
      return i;
  }
  return w.n_elem - 1;
}

static double row_log_sum_exp(const arma::rowvec &values)
{
  const double largest = values.max();
  if (!std::isfinite(largest))
    return largest;
  return largest + std::log(arma::accu(arma::exp(values - largest)));
}

double AveragedLatentStep::draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed,
                                      int direction, Particle* refreshed, bool* refreshed_made) const
{
  if (refreshed_made != NULL)
    *refreshed_made = false;
  // the variables the proposal changed (theta), and their proposed values
  std::vector<std::string> names;
  std::vector<arma::mat> end;
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
      if (a->first == this->variable)
        Rcpp::stop("averaged(…, latent = " + this->variable + "): the proposal changes " + this->variable + ", which the move takes along itself; propose the other parameters.");
      names.push_back(a->first);
      end.push_back(*a->second.first);
    }
  }
  // the current state at theta and at vartheta
  Particle at_theta = current;
  Particle at_vartheta = current;
  for (size_t j = 0; j < names.size(); ++j)
    at_vartheta.parameters.overwrite_variable(names[j], end[j]);
  remake_stores(at_vartheta);
  const double L_theta = current.get_cached_likelihoods(index);
  const double L_vartheta = at_vartheta.evaluate_likelihoods(index);
  if ((L_theta == -arma::datum::inf) || (L_vartheta == -arma::datum::inf) || std::isnan(L_vartheta))
    return -arma::datum::inf;
  // the candidates: z, and M - 1 draws from z's own distribution at q's parameters
  const Parameters &q_at = (direction == 1) ? at_vartheta.parameters : current.parameters;
  const arma::mat &z = static_cast<const Parameters&>(current.parameters)[this->variable];
  std::vector<arma::mat> candidates(this->points);
  candidates[0] = z;
  Parameters given = q_at;
  given.overwrite_variable(this->variable, arma::mat()); // (the draw writes z: not the particle's)
  for (size_t i = 1; i < this->points; ++i)
  {
    Parameters drawn = this->simulate(rng, given);
    const arma::mat* value = static_cast<const Parameters&>(drawn).find_vector(this->variable);
    if ((value == NULL) || (value->n_rows != z.n_rows) || (value->n_cols != z.n_cols))
      Rcpp::stop("averaged(…, latent = " + this->variable + "): a draw of " + this->variable + " from its distribution does not have its size.");
    candidates[i] = *value;
  }
  // the changes at each element and candidate, at theta and at vartheta
  arma::mat full_theta, own_theta, full_vartheta, own_vartheta;
  this->element_changes(at_theta, candidates, full_theta, own_theta);
  this->element_changes(at_vartheta, candidates, full_vartheta, own_vartheta);
  const arma::mat &D = (direction == 1) ? own_vartheta : own_theta;
  const arma::mat a = full_vartheta - D, b = full_theta - D;
  double log_ratio = L_vartheta - L_theta;
  for (size_t t = 0; t < a.n_rows; ++t)
    log_ratio += row_log_sum_exp(a.row(t)) - row_log_sum_exp(b.row(t));
  // the new latent, and the refreshed current one
  std::vector<size_t> chosen(a.n_rows), refreshed_choice(b.n_rows);
  for (size_t t = 0; t < a.n_rows; ++t)
    chosen[t] = draw_from(rng, a.row(t));
  if (refreshed != NULL)
    for (size_t t = 0; t < b.n_rows; ++t)
      refreshed_choice[t] = draw_from(rng, b.row(t));
  this->install(at_vartheta, candidates, chosen);
  proposed = std::move(at_vartheta);
  if (refreshed != NULL)
  {
    this->install(at_theta, candidates, refreshed_choice);
    *refreshed = std::move(at_theta);
    if (refreshed_made != NULL)
      *refreshed_made = true;
  }
  if (std::isnan(log_ratio))
    return -arma::datum::inf;
  return log_ratio;
}

}
