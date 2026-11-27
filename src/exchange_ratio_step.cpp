#include <algorithm>
#include "exchange_ratio_step.h"

namespace ilike
{

ExchangeRatioStep::ExchangeRatioStep(const std::vector<ExchangeField> &fields_in,
                                     Data* data_in,
                                     const std::string &tempering_variable_in)
{
  this->fields = fields_in;
  this->data = data_in;
  this->tempering_variable = tempering_variable_in;
}

RatioStep* ExchangeRatioStep::duplicate() const
{
  return new ExchangeRatioStep(*this);
}

double ExchangeRatioStep::draw(RandomNumberGenerator &rng,
                               const Parameters &from,
                               const Parameters &to) const
{
  double total = 0.0;
  // the fraction of the replicates in the target (data tempering)
  double fraction = 1.0;
  if (!this->tempering_variable.empty())
  {
    const std::vector<std::string> names = to.get_vector_variables();
    if (std::find(names.begin(), names.end(), this->tempering_variable) != names.end())
      fraction = to[this->tempering_variable][0];
  }
  for (const ExchangeField &field : this->fields)
  {
    double whole = fraction * double(field.replicates);
    size_t active = size_t(std::max(std::round(whole), 0.0));
    if (std::abs(whole - double(active)) > 1e-8)
      Rcpp::stop("exchange(…) under data tempering: the sequence must add whole replicates (data(…, anneal = false) with a schedule of "
                 "multiples of 1/replicates); the exchange estimate of a fraction of a normalising constant is not available.");
    for (size_t r = 0; r < active; ++r)
    {
    arma::mat u;
    if (field.engine)
    {
      // the general engine (a distribution written as a block): its kernel
      // at `to`, from the latent field or from a draw from the base measure
      Parameters state = to.deep_copy();
      if (state.find_vector(field.variable) == NULL)
        state[field.variable] = field.engine->draw_base(rng, to)[field.variable];
      field.engine->draw(rng, state, to, std::max(field.sweeps, size_t(1)));
      u = state[field.variable];
    }
    else if (field.exact != NULL)
    {
      Parameters drawn = field.exact(rng, to, *this->data);
      u = drawn[field.variable];
    }
    else
    {
      // sweeps (approximate): from the field in the parameters (a latent
      // field's current value), otherwise from the base measure
      Parameters state = to.deep_copy();
      for (size_t k = 0; k < std::max(field.sweeps, size_t(1)); ++k)
        field.sweep(rng, state, *this->data);
      u = state[field.variable];
    }
    // (a draw that could not be made, e.g. an exact ising draw at beta < 0:
    // density 0 there)
    if (u.has_nan())
      return -arma::datum::inf;
    // the exchange term as an annealed importance sampling estimate of
    // Z(from)/Z(to) along the path from `to` (where u was drawn) to `from`
    // (dev/plan.md 3.22 U3c): with no bridges, gamma(u | from)/gamma(u | to)
    // (Murray, Ghahramani & MacKay 2006); with bridges, their section 5, with
    // reversible updates (a sweep then a reverse sweep) at each bridge
    Parameters state = to.deep_copy();
    state[field.variable] = u;
    if (field.engine)
    {
      total += field.engine->single(rng, state, to, from);
      continue;
    }
    FieldPath path(field.variable, field.potential, field.sweep, field.reverse_sweep, field.normaliser,
                   field.bridges, field.bridge_sweeps, field.bridges > 0, this->data);
    total += path.single(rng, state, to, from);
    }
  }
  return total;
}

}
