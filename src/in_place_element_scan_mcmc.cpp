#include <cmath>
#include "in_place_element_scan_mcmc.h"
#include "esjd_adaptor.h"
#include "index.h"
#include "particle.h"
#include "vector_factor_variables.h"
#include "likelihood_estimator_output.h"
#include "plate_likelihood_estimator_output.h"
#include "exact_likelihood_estimator_output.h"
#include "annealed_likelihood_estimator_output.h"
#include "distributions.h"
#include "metropolis_mcmc.h"
#include "mcmc_adaptor.h"
#include "covariance_tools.h"
#include "smc_output.h"
#include "particles.h"
#include "move_output.h"
#include "ensemble_kalman_output.h"
#include "ensemble.h"

namespace ilike
{

InPlaceElementScanMCMC::InPlaceElementScanMCMC()
  :DeterministicScanMCMC()
{
  this->step_kind = ElementStep::normal;
}

InPlaceElementScanMCMC::InPlaceElementScanMCMC(MCMCTermination* termination_in,
                                               const std::vector<MCMC*> &element_moves_in,
                                               const std::string &variable_in,
                                               const std::vector<size_t> &rows_in,
                                               const std::vector<size_t> &columns_in,
                                               const std::vector<double> &scales_in,
                                               ElementStep step_in,
                                               const std::vector<double> &states_in)
  :DeterministicScanMCMC(termination_in, element_moves_in)
{
  this->variable = variable_in;
  this->rows = rows_in;
  this->columns = columns_in;
  this->scales = scales_in;
  this->step_kind = step_in;
  this->states = states_in;
}

InPlaceElementScanMCMC::~InPlaceElementScanMCMC()
{
}

InPlaceElementScanMCMC::InPlaceElementScanMCMC(const InPlaceElementScanMCMC &another)
  :DeterministicScanMCMC(another)
{
  this->make_copy(another);
}

void InPlaceElementScanMCMC::operator=(const InPlaceElementScanMCMC &another)
{
  if (this == &another)
    return;
  DeterministicScanMCMC::operator=(another);
  this->make_copy(another);
}

Kernel* InPlaceElementScanMCMC::duplicate() const
{
  return new InPlaceElementScanMCMC(*this);
}

MCMC* InPlaceElementScanMCMC::mcmc_duplicate() const
{
  return new InPlaceElementScanMCMC(*this);
}

DeterministicScanMCMC* InPlaceElementScanMCMC::deterministic_scan_mcmc_duplicate() const
{
  return new InPlaceElementScanMCMC(*this);
}

void InPlaceElementScanMCMC::make_copy(const InPlaceElementScanMCMC &another)
{
  this->variable = another.variable;
  this->rows = another.rows;
  this->columns = another.columns;
  this->scales = another.scales;
  this->step_kind = another.step_kind;
  this->states = another.states;
  this->smc_constant = another.smc_constant;
  this->index.reset((another.index != NULL) ? another.index->duplicate() : NULL);
}

ElementRandomWalkProposalKernel* InPlaceElementScanMCMC::element_kernel(size_t e) const
{
  if (e >= this->moves.size())
    return NULL;
  MetropolisMCMC* move = dynamic_cast<MetropolisMCMC*>(this->moves[e].get());
  if (move == NULL)
    return NULL;
  return dynamic_cast<ElementRandomWalkProposalKernel*>(move->get_symmetric_proposal());
}

// the element's step: its kernel's (which adaptation may change), or the
// scan's own for gibbs scans (which hold no element moves)
double InPlaceElementScanMCMC::step_scale(size_t e) const
{
  ElementRandomWalkProposalKernel* kernel = this->element_kernel(e);
  return (kernel != NULL) ? kernel->get_scale() : this->scales[e];
}

void InPlaceElementScanMCMC::adaptive_move_in_place(RandomNumberGenerator &rng,
                                                    Particle &particle,
                                                    Particle &scratch)
{
  // in place, adapting each element's kernel; otherwise the element moves'
  // own adaptive moves (the same kernels and adaptors)
  if (!this->scan_in_place(rng, particle, true))
    DeterministicScanMCMC::adaptive_move_in_place(rng, particle, scratch);
}

void InPlaceElementScanMCMC::append_adaptation_state(std::vector<std::string> *names,
                                                     std::vector<double> *values,
                                                     bool final,
                                                     const std::string &prefix) const
{
  // each element's columns name the element (theta[17].scale): no prefix per element move
  for (size_t k = 0; k < this->moves.size(); ++k)
    this->moves[k]->append_adaptation_state(names, values, final, prefix);
}

void InPlaceElementScanMCMC::set_steps_from(const arma::mat &points,
                                            const arma::colvec &normalised_weights,
                                            const std::vector<std::string> &names,
                                            const std::vector<size_t> &offsets,
                                            size_t number_of_rows)
{
  auto found = std::find(names.begin(), names.end(), this->variable);
  if ((found == names.end()) || (points.n_rows < 2))
    return;
  // each element's conditional variance given all the other parameters
  arma::colvec variances;
  if (!conditional_variances(points, normalised_weights, variances))
    return;
  size_t offset = offsets[found - names.begin()];
  for (size_t e = 0; e < this->rows.size(); ++e)
  {
    ElementRandomWalkProposalKernel* kernel = this->element_kernel(e);
    arma::uword column = offset + this->columns[e]*number_of_rows + this->rows[e];
    if ((kernel == NULL) || (column >= variances.n_elem))
      continue;
    double sd = std::sqrt(variances[column]);
    if (std::isfinite(sd) && (sd > 0.0))
      kernel->set_scale(this->smc_constant*sd);
  }
}

void InPlaceElementScanMCMC::smc_adapt(SMCOutput* current_state)
{
  if (this->smc_constant <= 0.0)
  {
    DeterministicScanMCMC::smc_adapt(current_state);
    return;
  }
  Particles &generation = current_state->back();
  std::vector<std::string> names;
  std::vector<size_t> offsets;
  arma::mat points = population_points(generation, names, offsets);
  if (points.n_rows == 0)
    return;
  arma::colvec weights;
  if (!population_weights(generation, points.n_rows, weights))
    return;
  size_t number_of_rows = generation.is_flat() ? generation.flat_points.parameters_at(0)[this->variable].n_rows
                                               : generation.particles[0]->back().parameters[this->variable].n_rows;
  this->set_steps_from(points, weights, names, offsets, number_of_rows);
}

void InPlaceElementScanMCMC::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  if (this->smc_constant <= 0.0)
  {
    DeterministicScanMCMC::ensemble_adapt(current_state);
    return;
  }
  Ensemble &ensemble = current_state->back();
  std::vector<std::string> names;
  std::vector<size_t> offsets;
  arma::mat points = ensemble_points(ensemble, names, offsets);
  if (points.n_rows == 0)
    return;
  size_t number_of_rows = ensemble.members[0]->back().parameters[this->variable].n_rows;
  this->set_steps_from(points, arma::colvec(points.n_rows, arma::fill::value(1.0/double(points.n_rows))),
                       names, offsets, number_of_rows);
}

void InPlaceElementScanMCMC::set_index(Index* index_in)
{
  DeterministicScanMCMC::set_index(index_in);
  this->index.reset((index_in != NULL) ? index_in->duplicate() : NULL);
  this->update_index_with_cut(this->index.get());
}

void InPlaceElementScanMCMC::set_index_if_null(Index* index_in)
{
  DeterministicScanMCMC::set_index_if_null(index_in);
  if ((this->index == NULL) && (index_in != NULL))
    this->index.reset(index_in->duplicate());
  this->update_index_with_cut(this->index.get());
}

void InPlaceElementScanMCMC::set_cut_index(const Index* cut_index_in)
{
  // (the element moves get theirs from the recipe)
  this->store_cut_index(cut_index_in);
  this->update_index_with_cut(this->index.get());
}

Particle InPlaceElementScanMCMC::move(RandomNumberGenerator &rng,
                                      const Particle &particle) const
{
  Particle moved = particle;
  this->move_in_place(rng, moved);
  return moved;
}

void InPlaceElementScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                           Particle &particle) const
{
  if (!this->scan_in_place(rng, particle))
  {
    if (this->is_gibbs())
      Rcpp::stop("gibbs moves need the model's factors held as a list of outputs (not available here).");
    DeterministicScanMCMC::move_in_place(rng, particle);
  }
}

void InPlaceElementScanMCMC::move_in_place(RandomNumberGenerator &rng,
                                           Particle &particle,
                                           Particle &scratch) const
{
  if (!this->scan_in_place(rng, particle))
  {
    if (this->is_gibbs())
      Rcpp::stop("gibbs moves need the model's factors held as a list of outputs (not available here).");
    DeterministicScanMCMC::move_in_place(rng, particle, scratch);
  }
}

bool InPlaceElementScanMCMC::scan_in_place(RandomNumberGenerator &rng,
                                           Particle &particle,
                                           bool adapt) const
{
  VectorFactorVariables* factor_variables = dynamic_cast<VectorFactorVariables*>(particle.factor_variables.get());
  if ((factor_variables == NULL) || (particle.ensemble_factor_variables != NULL) || (this->index == NULL))
    return false;

  // the variable's matrix owned by this particle alone (copied if shared, e.g.
  // with a stored iteration), so that it can be changed in place
  {
    Parameters own;
    own[this->variable] = static_cast<const Parameters&>(particle.parameters)[this->variable];
    particle.parameters.deep_overwrite_with_variables_in_argument(own);
  }
  arma::mat &x = particle.parameters[this->variable];

  const Index* move_index = this->index.get();
  const Index* new_outputs_index = this->get_new_outputs_index(move_index);
  std::vector<LikelihoodEstimatorOutput*> changed_plates; // (plates, or annealed plates, changed in place)
  std::vector<std::pair<size_t, std::unique_ptr<LikelihoodEstimatorOutput>>> replaced;
  bool any_accepted = false;
  size_t row = 0, column = 0;

  // the factors in `which`, given the element's new value: plates in place,
  // the others with new outputs (evaluated if `evaluate`; for gibbs they must
  // be exact)
  auto set_new_outputs = [&](const Index* which, bool evaluate) {
    changed_plates.clear();
    replaced.clear();
    for (auto i = which->begin(); i != which->end(); ++i)
    {
      LikelihoodEstimatorOutput* output = factor_variables->get_output(*i);
      if (output == NULL)
        continue;
      if (output->change_in_place(particle.parameters, this->variable, row, column))
      {
        changed_plates.push_back(output);
        continue;
      }
      if (this->is_gibbs() && (dynamic_cast<ExactLikelihoodEstimatorOutput*>(output) == NULL) &&
          (dynamic_cast<AnnealedLikelihoodEstimatorOutput*>(output) == NULL))
        Rcpp::stop("gibbs moves need exact factors: factor " + std::to_string(*i + 1) + " is estimated.");
      replaced.push_back(std::make_pair(*i, factor_variables->replace_with_new_output(*i, particle.parameters)));
      if (evaluate)
        factor_variables->get_output(*i)->evaluate(particle.parameters);
    }
  };
  // undo set_new_outputs
  auto restore = [&]() {
    for (auto p = changed_plates.rbegin(); p != changed_plates.rend(); ++p)
      (*p)->undo_in_place();
    for (auto r = replaced.rbegin(); r != replaced.rend(); ++r)
      factor_variables->put_back_output(r->first, std::move(r->second));
    changed_plates.clear();
    replaced.clear();
  };

  // the plates' terms valid at the current point, once for the scan (after
  // a restore from a summary in SMC they are not, and every element would
  // otherwise make and evaluate a new output)
  for (const Index* which : {move_index, new_outputs_index})
    for (auto i = which->begin(); i != which->end(); ++i)
    {
      LikelihoodEstimatorOutput* output = factor_variables->get_output(*i);
      if (output != NULL)
        output->prepare_in_place(particle.parameters);
    }

  // esjd (C226): the particle's row for the first element's kernel, which
  // holds the adaptor (one multiplier of every element's step)
  const ParticleTuning::Row* tuned_row = NULL;
  if ((particle.tuning != NULL) && !this->rows.empty())
  {
    ElementRandomWalkProposalKernel* first = this->element_kernel(0);
    if (first != NULL)
      for (const auto &r : particle.tuning->rows)
        if (r.instance == first->get_instance_index())
          tuned_row = &r;
  }

  for (size_t e_order = 0; e_order < this->rows.size(); ++e_order)
  {
    // (U5a) move_in_place_reversed: the elements in reverse order
    const size_t e = this->reversed_order ? (this->rows.size() - 1 - e_order) : e_order;
    row = this->rows[e];
    column = this->columns[e];
    if ((row >= x.n_rows) || (column >= x.n_cols))
      Rcpp::stop("The move on " + this->variable + "[" + std::to_string(row + 1) +
                 (x.n_cols > 1 || column > 0 ? ", " + std::to_string(column + 1) : std::string("")) +
                 "] is outside the variable, which has " + std::to_string(x.n_rows) + " row(s) and " +
                 std::to_string(x.n_cols) + " column(s).");
    double old_value = x(row, column);
    double denominator = particle.get_cached_likelihoods(move_index);

    // the proposal, and the factors' new outputs, in the order of the moves
    // through proposed particles (so the same random numbers)
    if (!this->states.empty())
    {
      // gibbs: the full conditional over the states
      std::vector<double> log_weights(this->states.size());
      double largest = -arma::datum::inf;
      for (size_t k = 0; k < this->states.size(); ++k)
      {
        x(row, column) = this->states[k];
        set_new_outputs(move_index, true);
        log_weights[k] = particle.get_cached_likelihoods(move_index);
        largest = std::max(largest, log_weights[k]);
        restore();
      }
      x(row, column) = old_value;
      double u = runif(rng);
      if (largest == -arma::datum::inf)
        continue;
      double total = 0.0;
      for (double w : log_weights)
        total = total + std::exp(w - largest);
      size_t chosen = this->states.size() - 1;
      double cumulative = 0.0;
      for (size_t k = 0; k < this->states.size(); ++k)
      {
        cumulative = cumulative + std::exp(log_weights[k] - largest)/total;
        if (u < cumulative)
        {
          chosen = k;
          break;
        }
      }
      x(row, column) = this->states[chosen];
      set_new_outputs(new_outputs_index, true);
      for (LikelihoodEstimatorOutput* plate : changed_plates)
        plate->commit_in_place();
      replaced.clear();
      particle.target_evaluated = log_weights[chosen];
      this->refresh_cut_factors(particle);
      this->check_carried_over_factors(particle, move_index);
      if (this->states[chosen] != old_value)
        any_accepted = true;
      continue;
    }

    if (this->step_kind == ElementStep::normal_conditional)
    {
      // the log target as a function of the element, f(v) = A v^2 + B v + C
      // when the conditional is normal: from f at old - h, old, old + h; h is
      // then set to the conditional's sd and f evaluated again if they differ
      // by more than a factor of 10 (the second difference loses relative
      // accuracy about eps |f| / (|A| h^2))
      auto target_at = [&](double v) {
        x(row, column) = v;
        set_new_outputs(move_index, true);
        double value = particle.get_cached_likelihoods(move_index);
        restore();
        return value;
      };
      double centre = denominator;
      double h = std::max(1.0, std::abs(old_value))*1e-2;
      double mean = 0.0, sd = 0.0;
      for (int pass = 0; pass < 2; ++pass)
      {
        double below = target_at(old_value - h), above = target_at(old_value + h);
        double A = (above - 2.0*centre + below)/(2.0*h*h);
        double slope = (above - below)/(2.0*h); // f'(old) = 2 A old + B
        if (!(A < 0.0) || !std::isfinite(A) || !std::isfinite(slope))
          Rcpp::stop("gibbs(" + this->variable + "[…]): the log density of " + this->variable + "[" +
                     std::to_string(row + 1) + (x.n_cols > 1 ? ", " + std::to_string(column + 1) : std::string("")) +
                     "] given the rest is not a concave quadratic at its current value, so its full conditional is not normal "
                     "(is every factor reading it normal in it?).");
        sd = std::sqrt(-1.0/(2.0*A));
        mean = old_value - slope/(2.0*A);
        if ((sd < 10.0*h) && (sd > 0.1*h))
          break;
        h = sd;
      }
      x(row, column) = old_value;
      double proposed = rnorm(rng, mean, sd);
      x(row, column) = proposed;
      set_new_outputs(new_outputs_index, false);
      double log_u = log(runif(rng));
      for (auto &r : replaced)
        if (std::find(move_index->begin(), move_index->end(), r.first) != move_index->end())
          factor_variables->get_output(r.first)->evaluate(particle.parameters);
      double numerator = particle.get_cached_likelihoods(move_index);
      // independent proposal q = N(mean, sd^2): ratio pi(new) q(old) / (pi(old) q(new))
      double log_q_old = -0.5*(old_value - mean)*(old_value - mean)/(sd*sd);
      double log_q_new = -0.5*(proposed - mean)*(proposed - mean)/(sd*sd);
      bool accept = (numerator != -arma::datum::inf) && (denominator != -arma::datum::inf) &&
                    (log_u < numerator - denominator + log_q_old - log_q_new);
      if (accept)
      {
        for (LikelihoodEstimatorOutput* plate : changed_plates)
          plate->commit_in_place();
        replaced.clear();
        particle.target_evaluated = numerator;
        this->refresh_cut_factors(particle);
        this->check_carried_over_factors(particle, move_index);
        any_accepted = true;
      }
      else
      {
        x(row, column) = old_value;
        restore();
      }
      continue;
    }

    if ((this->step_kind == ElementStep::gamma_conditional) || (this->step_kind == ElementStep::invgamma_conditional) ||
        (this->step_kind == ElementStep::beta_conditional))
    {
      // the log target as a function of the element, f(v) = a phi1(v) +
      // b phi2(v) + c exactly when the conditional is in the family (C215):
      // from f at three values inside the support (around the current value
      // on the log or logit scale), a and b by a 2 x 2 solve
      bool beta = (this->step_kind == ElementStep::beta_conditional);
      bool inverse = (this->step_kind == ElementStep::invgamma_conditional);
      auto phi = [&](double v, double &p1, double &p2) {
        p1 = std::log(v);
        p2 = beta ? std::log1p(-v) : (inverse ? 1.0/v : v);
      };
      auto target_at = [&](double v) {
        x(row, column) = v;
        set_new_outputs(move_index, true);
        double value = particle.get_cached_likelihoods(move_index);
        restore();
        return value;
      };
      bool inside = beta ? ((old_value > 0.0) && (old_value < 1.0)) : (old_value > 0.0);
      double v[3], f[3];
      const double h = 0.5;
      if (beta)
      {
        double centre = inside ? std::log(old_value) - std::log1p(-old_value) : 0.0;
        for (int k = 0; k < 3; ++k)
          v[k] = 1.0/(1.0 + std::exp(-(centre + double(k - 1)*h)));
      }
      else
      {
        double centre = inside ? old_value : 1.0;
        for (int k = 0; k < 3; ++k)
          v[k] = centre*std::exp(double(k - 1)*h);
      }
      for (int k = 0; k < 3; ++k)
        f[k] = ((k == 1) && inside) ? denominator : target_at(v[k]);
      double p[3][2];
      for (int k = 0; k < 3; ++k)
        phi(v[k], p[k][0], p[k][1]);
      double m00 = p[0][0] - p[1][0], m01 = p[0][1] - p[1][1], m10 = p[2][0] - p[1][0], m11 = p[2][1] - p[1][1];
      double r0 = f[0] - f[1], r1 = f[2] - f[1];
      double determinant = m00*m11 - m01*m10;
      double a = (r0*m11 - m01*r1)/determinant, b = (m00*r1 - r0*m10)/determinant;
      // the conditional's parameters
      double first, second;
      const char* family = beta ? "beta" : (inverse ? "inverse gamma" : "gamma");
      if (beta)
      {
        first = a + 1.0;
        second = b + 1.0;
      }
      else if (inverse)
      {
        first = -a - 1.0; // shape
        second = -b;      // scale
      }
      else
      {
        first = a + 1.0;  // shape
        second = -b;      // rate
      }
      if (!(first > 0.0) || !(second > 0.0) || !std::isfinite(first) || !std::isfinite(second))
        Rcpp::stop("gibbs(" + this->variable + "[…]): the log density of " + this->variable + "[" +
                   std::to_string(row + 1) + (x.n_cols > 1 ? ", " + std::to_string(column + 1) : std::string("")) +
                   "] given the rest is not that of a " + family + " distribution (fitted parameters " +
                   std::to_string(first) + ", " + std::to_string(second) + ").");
      x(row, column) = old_value;
      double proposed;
      if (beta)
      {
        double g1 = rgamma(rng, first, 1.0), g2 = rgamma(rng, second, 1.0);
        proposed = g1/(g1 + g2);
      }
      else if (inverse)
        proposed = 1.0/rgamma(rng, first, second);
      else
        proposed = rgamma(rng, first, second);
      x(row, column) = proposed;
      set_new_outputs(new_outputs_index, false);
      double log_u = log(runif(rng));
      for (auto &r : replaced)
        if (std::find(move_index->begin(), move_index->end(), r.first) != move_index->end())
          factor_variables->get_output(r.first)->evaluate(particle.parameters);
      double numerator = particle.get_cached_likelihoods(move_index);
      // independent proposal q(v) proportional to exp(a phi1(v) + b phi2(v))
      // (the density's own log-Jacobian is inside a phi1)
      double q_old = -arma::datum::inf, q_new, p1, p2;
      if (inside)
      {
        phi(old_value, p1, p2);
        q_old = a*p1 + b*p2;
      }
      phi(proposed, p1, p2);
      q_new = a*p1 + b*p2;
      // (from a current value outside the support, e.g. a starting value, the
      // draw is accepted: the target there is 0)
      bool accept = (numerator != -arma::datum::inf) && std::isfinite(q_new) &&
                    ((denominator == -arma::datum::inf) || (log_u < numerator - denominator + q_old - q_new));
      if (accept)
      {
        for (LikelihoodEstimatorOutput* plate : changed_plates)
          plate->commit_in_place();
        replaced.clear();
        particle.target_evaluated = numerator;
        this->refresh_cut_factors(particle);
        this->check_carried_over_factors(particle, move_index);
        any_accepted = true;
      }
      else
      {
        x(row, column) = old_value;
        restore();
      }
      continue;
    }

    ElementRandomWalkProposalKernel* kernel = adapt ? this->element_kernel(e) : NULL;
    if (kernel != NULL)
      kernel->prepare_mcmc_adaptation();
    double scale = this->step_scale(e);
    // a particle's own multiplier of the steps (esjd in SMC, C226)
    if (tuned_row != NULL)
      scale = scale*tuned_row->values[0];
    double proposed = (this->step_kind == ElementStep::normal) ? rnorm(rng, old_value, scale)
                    : (this->step_kind == ElementStep::uniform) ? runif(rng, old_value - scale, old_value + scale)
                    : 1.0 - old_value;
    x(row, column) = proposed;
    set_new_outputs(new_outputs_index, false);
    double log_u = log(runif(rng));

    // the numerator: the new outputs in the move's index evaluated (the
    // plates already hold their new values)
    for (auto &r : replaced)
      if (std::find(move_index->begin(), move_index->end(), r.first) != move_index->end())
        factor_variables->get_output(r.first)->evaluate(particle.parameters);
    double numerator = particle.get_cached_likelihoods(move_index);

    bool accept = (numerator != -arma::datum::inf) && (denominator != -arma::datum::inf) &&
                  (log_u < numerator - denominator);
    if (tuned_row != NULL)
    {
      // the expected squared jump of this step, alpha (y - x)^2 M_kk
      double alpha = 0.0;
      if ((numerator != -arma::datum::inf) && (denominator != -arma::datum::inf))
      {
        double log_ratio = numerator - denominator;
        alpha = (log_ratio >= 0.0) ? 1.0 : (std::isnan(log_ratio) ? 0.0 : std::exp(log_ratio));
      }
      double delta = proposed - old_value;
      double jump = alpha*delta*delta*tuned_row->adaptor->metric_at(row + column*x.n_rows);
      tuned_row->jumps = tuned_row->jumps + (std::isfinite(jump) ? jump : 0.0);
      tuned_row->iterations = tuned_row->iterations + 1;
    }
    if (kernel != NULL)
    {
      // the adaptation step (as MetropolisMCMC::adaptive_move_in_place):
      // alpha = min(1, ratio), 0 when either density is 0; the element
      // adaptors read only alpha
      double acceptance_probability = 0.0;
      if ((numerator != -arma::datum::inf) && (denominator != -arma::datum::inf))
      {
        double log_ratio = numerator - denominator;
        if (log_ratio >= 0.0)
          acceptance_probability = 1.0;
        else if (!std::isnan(log_ratio))
          acceptance_probability = std::exp(log_ratio);
      }
      AdaptationStep step;
      step.previous = &particle;
      step.proposed = &particle;
      step.current = &particle;
      step.acceptance_probability = acceptance_probability;
      kernel->mcmc_adapt(step);
    }
    if (accept)
    {
      for (LikelihoodEstimatorOutput* plate : changed_plates)
        plate->commit_in_place();
      replaced.clear();
      particle.target_evaluated = numerator;
      this->refresh_cut_factors(particle);
      this->check_carried_over_factors(particle, move_index);
      any_accepted = true;
    }
    else
    {
      x(row, column) = old_value;
      restore();
    }
  }

  // what the particle keeps for proposals (transformed parameters, gradients)
  // was for its old parameters: made again, as for a proposed particle
  if (any_accepted)
    particle.simulate_proposal_variables(particle.proposals_to_transform_for_pointer,
                                         particle.proposals_to_find_gradient_for_pointer);
  return true;
}


void InPlaceElementScanMCMC::move_in_place_reversed(RandomNumberGenerator &rng, Particle &particle) const
{
  this->reversed_order = true;
  bool done = this->scan_in_place(rng, particle);
  this->reversed_order = false;
  if (!done)
  {
    if (this->is_gibbs())
      Rcpp::stop("gibbs moves need the model's factors held as a list of outputs (not available here).");
    DeterministicScanMCMC::move_in_place_reversed(rng, particle);
  }
}

}
