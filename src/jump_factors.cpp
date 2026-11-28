#include <algorithm>
#include "jump_factors.h"
#include "likelihood_estimator.h"
#include "likelihood_estimator_output.h"

namespace ilike
{

const char* jump_prefix = "ilike_jump_";

static bool prefixed(const std::string &name)
{
  return name.compare(0, std::string(jump_prefix).size(), jump_prefix) == 0;
}

Parameters jump_model_parameters(const Parameters &input)
{
  Parameters output;
  for (auto v = input.vector_begin(); v != input.vector_end(); ++v)
    if (!prefixed(v->first))
      output[v->first] = *v->second.first;
  return output;
}

Parameters jump_reverse_parameters(const Parameters &input)
{
  Parameters output;
  const size_t n = std::string(jump_prefix).size();
  for (auto v = input.vector_begin(); v != input.vector_end(); ++v)
    output[prefixed(v->first) ? v->first.substr(n) : v->first] = *v->second.first;
  return output;
}

double jump_start_power(const Parameters &parameters, const std::string &variable)
{
  const arma::mat* t = parameters.find_vector(variable);
  return (t != NULL) ? 1.0 - (*t)(0) : 0.0;
}

double jump_end_power(const Parameters &parameters, const std::string &variable)
{
  const arma::mat* t = parameters.find_vector(variable);
  return (t != NULL) ? (*t)(0) : 0.0;
}

std::string jump_path_variable(size_t jump)
{
  return "ilike_jump_t_" + std::to_string(jump);
}

// ---- q(u' | x') ---------------------------------------------------------------

ReverseAuxiliaryFactor::ReverseAuxiliaryFactor(InvolutionFunctions* functions_in, Data* data_in)
  : LikelihoodFactor(data_in), functions(functions_in)
{
}

ReverseAuxiliaryFactor::ReverseAuxiliaryFactor(const ReverseAuxiliaryFactor &another)
  : LikelihoodFactor(another), functions(another.functions->duplicate())
{
}

ReverseAuxiliaryFactor::~ReverseAuxiliaryFactor()
{
}

Factor* ReverseAuxiliaryFactor::duplicate() const
{
  return new ReverseAuxiliaryFactor(*this);
}

LikelihoodFactor* ReverseAuxiliaryFactor::likelihood_factor_duplicate() const
{
  return new ReverseAuxiliaryFactor(*this);
}

double ReverseAuxiliaryFactor::likelihood_evaluate(const Parameters &input) const
{
  return this->functions->log_density(jump_reverse_parameters(input));
}

arma::mat ReverseAuxiliaryFactor::likelihood_evaluate_gradient(const std::string &variable, const Parameters &input) const
{
  Rcpp::stop("Annealed reversible jump: gradients of the bridges are not available yet, so the kernel's moves cannot use gradients.");
}

// ---- f0 -------------------------------------------------------------------------

MappedBackFactor::MappedBackFactor(InvolutionFunctions* functions_in, const std::string &indicator_in,
                                   const std::vector<std::string> &variables_in,
                                   const std::vector<LikelihoodEstimator*> &inner_in, Data* data_in)
  : LikelihoodFactor(data_in), functions(functions_in), indicator(indicator_in), variables(variables_in)
{
  for (LikelihoodEstimator* e : inner_in)
    this->inner.emplace_back(e);
}

MappedBackFactor::MappedBackFactor(const MappedBackFactor &another)
  : LikelihoodFactor(another), functions(another.functions->duplicate()), indicator(another.indicator),
    variables(another.variables)
{
  for (const auto &e : another.inner)
    this->inner.emplace_back(e->duplicate());
}

MappedBackFactor::~MappedBackFactor()
{
}

Factor* MappedBackFactor::duplicate() const
{
  return new MappedBackFactor(*this);
}

LikelihoodFactor* MappedBackFactor::likelihood_factor_duplicate() const
{
  return new MappedBackFactor(*this);
}

double MappedBackFactor::likelihood_evaluate(const Parameters &input) const
{
  const std::string to_name = this->indicator + "_proposed";
  // xi' = (k', z', k, u'), and (z, u) = phi(xi')
  Parameters reverse = jump_reverse_parameters(input);
  Parameters image = this->functions->map(reverse);
  // the state at k
  Parameters mapped = jump_model_parameters(input);
  const double k_new = static_cast<const Parameters&>(mapped)[this->indicator](0);
  for (const std::string &v : this->variables)
  {
    const arma::mat* value = static_cast<const Parameters&>(image).find_vector(v);
    if (value == NULL)
      Rcpp::stop("Annealed reversible jump: the map does not return " + v + ".");
    mapped[v] = *value;
  }
  mapped[this->indicator] = static_cast<const Parameters&>(reverse)[to_name];
  // xi = (k, z, k', u): u, what the map returns besides the jump's variables
  Parameters source = mapped;
  source[to_name] = arma::mat(1, 1, arma::fill::value(k_new));
  for (auto v = image.vector_begin(); v != image.vector_end(); ++v)
    if ((std::find(this->variables.begin(), this->variables.end(), v->first) == this->variables.end()) &&
        (v->first != this->indicator) && (v->first != to_name))
      source[v->first] = *v->second.first;
  // pi(k, z): the move's factors
  if (this->inner_outputs.size() != this->inner.size())
  {
    this->inner_outputs.clear();
    for (const auto &e : this->inner)
    {
      e->setup(mapped);
      this->inner_outputs.emplace_back(e->initialise(mapped));
      this->inner_outputs.back()->simulate(mapped);
    }
  }
  double log_pi = 0.0;
  for (const auto &output : this->inner_outputs)
  {
    log_pi += output->evaluate(mapped);
    if (log_pi == -arma::datum::inf)
      return log_pi;
  }
  return log_pi + this->functions->log_density(source) - this->functions->log_jacobian(source);
}

arma::mat MappedBackFactor::likelihood_evaluate_gradient(const std::string &variable, const Parameters &input) const
{
  Rcpp::stop("Annealed reversible jump: gradients of the bridges are not available yet, so the kernel's moves cannot use gradients.");
}

}
