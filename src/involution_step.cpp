#include <cmath>
#include <algorithm>
#include "involution_step.h"
#include "particle.h"
#include "index.h"
#include "distributions.h"
#include "factor_variables.h"
#include "jump_factors.h"

namespace ilike
{

// ---- the user's functions (involution, rj) -----------------------------------

UserInvolutionFunctions::UserInvolutionFunctions(InvolutionAuxiliaryPtr auxiliary_in, InvolutionLogDensityPtr log_auxiliary_in,
                                                 InvolutionMapPtr map_in, InvolutionLogDensityPtr log_jacobian_in, const Data* data_in)
{
  this->auxiliary = auxiliary_in;
  this->log_auxiliary = log_auxiliary_in;
  this->map_function = map_in;
  this->log_jacobian_function = log_jacobian_in;
  this->data = data_in;
}

Parameters UserInvolutionFunctions::draw(RandomNumberGenerator &rng, const Parameters &xi) const
{
  return this->auxiliary(rng, xi, *this->data);
}

double UserInvolutionFunctions::log_density(const Parameters &xi) const
{
  return this->log_auxiliary(xi, *this->data);
}

Parameters UserInvolutionFunctions::map(const Parameters &xi) const
{
  return this->map_function(xi, *this->data);
}

double UserInvolutionFunctions::log_jacobian(const Parameters &xi) const
{
  return (this->log_jacobian_function != NULL) ? this->log_jacobian_function(xi, *this->data) : 0.0;
}

// ---- birth_death's functions -------------------------------------------------

static const std::string birth_position = "ilike_birth_position";
static const std::string death_position = "ilike_death_position";

static std::string born_name(const std::string &variable)
{
  return "ilike_birth_" + variable;
}

BirthDeathFunctions::BirthDeathFunctions(const std::string &indicator_in, size_t maximum_in, const std::vector<std::string> &variables_in,
                                         SimulateIndependentProposalPtr simulate_in, EvaluateLogIndependentProposalPtr log_density_in)
{
  this->indicator = indicator_in;
  this->maximum = maximum_in;
  this->variables = variables_in;
  this->simulate = simulate_in;
  this->log_density_function = log_density_in;
}

static double scalar(const Parameters &xi, const std::string &name)
{
  return static_cast<const Parameters&>(xi)[name](0);
}

Parameters BirthDeathFunctions::draw(RandomNumberGenerator &rng, const Parameters &xi) const
{
  const double k = scalar(xi, this->indicator), k_new = scalar(xi, this->indicator + "_proposed");
  const size_t K = size_t(std::lround(k));
  Parameters u;
  if (k_new > k)
  {
    // a birth at position j (0-based, 0..K), and the new entries
    const size_t j = std::min(size_t(runif(rng) * double(K + 1)), K);
    u[birth_position] = arma::mat(1, 1, arma::fill::value(double(j)));
    Parameters born = this->simulate(rng, xi);
    for (const std::string &v : this->variables)
    {
      const arma::mat* b = static_cast<const Parameters&>(born).find_vector(v);
      if ((b == NULL) || (b->n_elem != 1))
        Rcpp::stop("birth_death: the new entry of " + v + " was not drawn as one value.");
      u[born_name(v)] = *b;
    }
  }
  else
  {
    // a death of entry j (0-based, 0..K-1)
    const size_t j = std::min(size_t(runif(rng) * double(K)), K - 1);
    u[death_position] = arma::mat(1, 1, arma::fill::value(double(j)));
  }
  return u;
}

double BirthDeathFunctions::log_density(const Parameters &xi) const
{
  // (the position's probability, 1/(k + 1) for a birth from k and for the
  // death back from k + 1, cancels: left out on both sides)
  if (static_cast<const Parameters&>(xi).find_vector(birth_position) == NULL)
    return 0.0;
  Parameters born;
  for (const std::string &v : this->variables)
    born[v] = static_cast<const Parameters&>(xi)[born_name(v)];
  return this->log_density_function(born, xi);
}

Parameters BirthDeathFunctions::map(const Parameters &xi) const
{
  const double k = scalar(xi, this->indicator);
  const size_t K = size_t(std::lround(k));
  const bool birth = (static_cast<const Parameters&>(xi).find_vector(birth_position) != NULL);
  const size_t j = size_t(std::lround(scalar(xi, birth ? birth_position : death_position)));
  Parameters image;
  for (const std::string &v : this->variables)
  {
    // (padded to the indicator's largest value)
    arma::mat values = arma::vectorise(static_cast<const Parameters&>(xi)[v]);
    if (values.n_elem < this->maximum)
      values.resize(this->maximum, 1);
    if (birth)
    {
      for (size_t i = K; i > j; --i)
        values(i) = values(i - 1);
      values(j) = scalar(xi, born_name(v));
    }
    else
    {
      image[born_name(v)] = arma::mat(1, 1, arma::fill::value(values(j)));
      for (size_t i = j; i + 1 < K; ++i)
        values(i) = values(i + 1);
    }
    image[v] = values;
  }
  // u': the reverse move's position
  image[birth ? death_position : birth_position] = arma::mat(1, 1, arma::fill::value(double(j)));
  return image;
}

bool BirthDeathFunctions::possible(double k, double k_new) const
{
  return ((k_new == k + 1.0) || (k_new == k - 1.0)) && (k_new >= 1.0) && (k_new <= double(this->maximum));
}

// ---- the step ----------------------------------------------------------------

InvolutionStep::InvolutionStep(InvolutionFunctions* functions_in, const std::string &indicator_in, const std::vector<std::string> &variables_in, size_t check_in,
                               const AnnealedPathKernel &annealing_in)
  : annealing(annealing_in)
{
  this->functions.reset(functions_in);
  this->indicator = indicator_in;
  this->variables = variables_in;
  this->check = check_in;
  this->checked_draws = std::make_shared<std::atomic<size_t>>(0);
}

InvolutionStep::InvolutionStep(const InvolutionStep &another)
  : annealing(another.annealing)
{
  this->functions.reset(another.functions->duplicate());
  this->indicator = another.indicator;
  this->variables = another.variables;
  this->check = another.check;
  this->checked_draws = another.checked_draws; // (one count for the move and its copies)
}

InvolutionStep::~InvolutionStep()
{
}

RatioStep* InvolutionStep::duplicate() const
{
  return new InvolutionStep(*this);
}

double InvolutionStep::draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const
{
  Rcpp::stop("InvolutionStep::draw: an involution move's step changes the state (draw_state).");
}

double InvolutionStep::first_direction_probability(const Particle &from, const Particle &to) const
{
  return this->functions->first_direction_probability(scalar(from.parameters, this->indicator), scalar(to.parameters, this->indicator));
}

static double largest_difference(const arma::mat &a, const arma::mat &b)
{
  if (arma::size(a) != arma::size(b))
    return arma::datum::inf;
  if (a.n_elem == 0)
    return 0.0;
  return arma::abs(a - b).max() / (1.0 + std::max(arma::abs(a).max(), arma::abs(b).max()));
}

static bool contains(const std::vector<std::string> &names, const std::string &name)
{
  return std::find(names.begin(), names.end(), name) != names.end();
}

double InvolutionStep::draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const
{
  const bool jumps = !this->indicator.empty();
  const std::string move_name = jumps ? "rj(…)" : "involution(…)";
  const std::string to_name = this->indicator + "_proposed";
  double k = 0.0, k_new = 0.0;
  if (jumps)
  {
    k = scalar(current.parameters, this->indicator);
    k_new = scalar(proposed.parameters, this->indicator);
    if (!this->functions->possible(k, k_new))
      return -arma::datum::inf;
  }
  // xi = (x, k', u): the state, the part of the auxiliary variables the
  // move's proposal chose (k'), and the rest, u ~ q(. | x, k')
  Parameters xi = current.parameters.deep_copy();
  if (jumps)
    xi[to_name] = arma::mat(1, 1, arma::fill::value(k_new));
  Parameters u = this->functions->draw(rng, xi);
  const std::vector<std::string> auxiliary_names = u.get_vector_variables();
  const std::vector<std::string> existing = current.parameters.get_vector_variables();
  for (const std::string &name : auxiliary_names)
  {
    if (contains(existing, name))
      Rcpp::stop(move_name + ": the auxiliary variable " + name + " has the name of a variable of the model.");
    xi[name] = static_cast<const Parameters&>(u)[name];
  }
  const double log_q = this->functions->log_density(xi);
  const double log_jacobian = this->functions->log_jacobian(xi);
  // (x', u') = phi(x, u); the new state is the proposed particle (its
  // factors have new outputs) with the variables the map changes
  Parameters image = this->functions->map(xi);
  Particle out = proposed;
  for (const std::string &name : this->variables)
  {
    const arma::mat* value = static_cast<const Parameters&>(image).find_vector(name);
    if (value == NULL)
      Rcpp::stop(move_name + ": the map does not return " + name + ".");
    out.parameters.overwrite_variable(name, *value);
  }
  Parameters xi_image = out.parameters.deep_copy();
  if (jumps)
    xi_image[to_name] = arma::mat(1, 1, arma::fill::value(k));
  // u': what the map returns besides the model's variables (it may differ
  // from u: a death draws only a position, and gives the removed values)
  for (const std::string &name : image.get_vector_variables())
  {
    if (contains(existing, name) || (name == to_name))
      continue;
    xi_image[name] = static_cast<const Parameters&>(image)[name];
  }
  // the check that the map undoes itself, on the first draws
  if (this->functions->checked() && (this->checked_draws->load() < this->check))
  {
    (*this->checked_draws)++;
    const std::string wrong = jumps ? "rj(…): the map for the jump back does not undo the map for the jump: applying both does not give back "
                                    : "involution(…): the map is not an involution: applying it twice does not give back ";
    Parameters back = this->functions->map(xi_image);
    for (const std::string &name : this->variables)
    {
      const arma::mat &before = static_cast<const Parameters&>(xi)[name];
      const arma::mat* again = static_cast<const Parameters&>(back).find_vector(name);
      bool same;
      if (jumps)
      {
        // (only the entries in the model before the jump)
        const size_t active = std::min(size_t(std::lround(k)), size_t(before.n_elem));
        same = (again != NULL) && (again->n_elem >= active) &&
          (largest_difference(arma::mat(arma::vectorise(*again)).head_rows(active), arma::mat(arma::vectorise(before)).head_rows(active)) <= 1e-8);
      }
      else
        same = (again != NULL) && (largest_difference(*again, before) <= 1e-8);
      if (!same)
        Rcpp::stop(wrong + name + ".");
    }
    for (const std::string &name : auxiliary_names)
    {
      const arma::mat* again = static_cast<const Parameters&>(back).find_vector(name);
      if ((again == NULL) || (largest_difference(*again, static_cast<const Parameters&>(xi)[name]) > 1e-8))
        Rcpp::stop(wrong + "the auxiliary variable " + name + ".");
    }
  }
  const double log_q_reverse = this->functions->log_density(xi_image);
  out.current_proposal_store.clear();
  out.simulate_proposal_variables(out.proposals_to_transform_for_pointer, out.proposals_to_find_gradient_for_pointer);
  // (the proposed particle's outputs carried over from the current ones what
  // the proposal left unchanged; the map changed more)
  if ((out.factor_variables != NULL) && (current.factor_variables != NULL))
    out.factor_variables->redo_carry_over(current.factor_variables.get(), out.parameters, index);
  const double L = current.get_cached_likelihoods(index);
  const double L_new = out.evaluate_likelihoods(index);
  if ((L_new == -arma::datum::inf) || std::isnan(L_new) || !std::isfinite(log_q) || std::isnan(log_q_reverse))
  {
    proposed = std::move(out);
    return -arma::datum::inf;
  }
  if (!this->annealing.runs())
  {
    proposed = std::move(out);
    // the target part: log pi(x') - log pi(x) + log q(u' | x') - log q(u | x) + log |J|
    return L_new - L + log_q_reverse - log_q + log_jacobian;
  }
  // (stage B) annealed reversible jump: the path from f0 (the model at k
  // through the map) to f1 (the model at k'), starting at x' with the reverse
  // move's auxiliary variables (and k, as <k>_proposed) under prefixed names;
  // its weight is the target part (with no bridges it is the one above)
  this->annealing.refuse_estimates(current, index, jumps ? "rj" : "involution");
  // the path runs in the larger model's coordinates, with the kernel moving
  // that model's variables and the auxiliary variables of that side fixed:
  // for a jump up from the new state x' = phi(x, u) (with u' and k), for a
  // jump down from the current state x (with u and k'), its end then mapped.
  // So the reverse of a jump's path is the reverse jump's path read backwards
  const bool up = (k_new > k);
  const Particle &side = up ? out : current;
  const Parameters &side_xi = up ? xi_image : xi;
  Parameters path_state = side.parameters;
  for (auto v = side_xi.vector_begin(); v != side_xi.vector_end(); ++v)
    if (side.parameters.find_vector(v->first) == NULL)
      path_state[std::string(jump_prefix) + v->first] = *v->second.first;
  const SMCPath* engine = up ? this->annealing.path.get() : this->annealing.reverse_path.get();
  if (engine == NULL)
    Rcpp::stop("rj(…, bridges = K): the path's SMC was not built.");
  const double log_weight = engine->single_state(rng, path_state, std::vector<std::string>(),
                                                 Parameters(), Parameters(), current.tuning);
  if (std::isnan(log_weight) || (log_weight == -arma::datum::inf))
  {
    proposed = std::move(out);
    return -arma::datum::inf;
  }
  if (up)
  {
    for (const std::string &name : this->variables)
      out.parameters.overwrite_variable(name, static_cast<const Parameters&>(path_state)[name]);
  }
  else
  {
    // the path's end, mapped to the smaller model
    Parameters end_image = this->functions->map(jump_reverse_parameters(path_state));
    for (const std::string &name : this->variables)
    {
      const arma::mat* value = static_cast<const Parameters&>(end_image).find_vector(name);
      if (value == NULL)
        Rcpp::stop(move_name + ": the map does not return " + name + ".");
      out.parameters.overwrite_variable(name, *value);
    }
  }
  // (the path changed the jump's variables: every factor at the new state)
  if ((out.factor_variables != NULL) && (current.factor_variables != NULL))
    out.factor_variables->redo_carry_over(current.factor_variables.get(), out.parameters, this->annealing.state_index.get());
  out.current_proposal_store.clear();
  out.simulate_proposal_variables(out.proposals_to_transform_for_pointer, out.proposals_to_find_gradient_for_pointer);
  const double all = out.evaluate_likelihoods(this->annealing.state_index.get());
  proposed = std::move(out);
  if ((all == -arma::datum::inf) || std::isnan(all))
    return -arma::datum::inf;
  return log_weight;
}

}
