#include "elementwise_transform_component.h"

namespace ilike
{

ElementwiseTransformComponent::ElementwiseTransformComponent()
:TransformComponent()
{
}

ElementwiseTransformComponent::ElementwiseTransformComponent(const std::vector<Spec> &specs_in)
:TransformComponent(), specs(specs_in)
{
  for (const Spec &spec : this->specs)
    if ((spec.kind == Kind::logit) && !(spec.lower < spec.upper))
      Rcpp::stop("logit(" + spec.variable + ", a, b): a must be below b.");
}

ElementwiseTransformComponent::~ElementwiseTransformComponent()
{
}

ElementwiseTransformComponent::ElementwiseTransformComponent(const ElementwiseTransformComponent &another)
:TransformComponent(another), specs(another.specs)
{
}

TransformComponent* ElementwiseTransformComponent::duplicate() const
{
  return new ElementwiseTransformComponent(*this);
}

const ElementwiseTransformComponent::Spec* ElementwiseTransformComponent::find(const std::string &variable) const
{
  for (const Spec &spec : this->specs)
    if (spec.variable == variable)
      return &spec;
  return NULL;
}

double ElementwiseTransformComponent::forward(const Spec &spec, double x)
{
  switch (spec.kind)
  {
    case Kind::log: return std::log(x);
    case Kind::lower: return std::log(x - spec.lower);
    case Kind::upper: return std::log(spec.upper - x);
    case Kind::logit:
    {
      double u = (x - spec.lower)/(spec.upper - spec.lower);
      return std::log(u) - std::log1p(-u);
    }
  }
  return x;
}

double ElementwiseTransformComponent::backward(const Spec &spec, double z)
{
  switch (spec.kind)
  {
    case Kind::log: return std::exp(z);
    case Kind::lower: return spec.lower + std::exp(z);
    case Kind::upper: return spec.upper - std::exp(z);
    case Kind::logit: return spec.lower + (spec.upper - spec.lower)/(1.0 + std::exp(-z));
  }
  return z;
}

double ElementwiseTransformComponent::log_dx_dz(const Spec &spec, double z)
{
  if (spec.kind != Kind::logit)
    return z; // (|dx/dz| = e^z)
  // log((b - a) s (1 - s)), s = 1/(1 + e^-z), computed stably
  return std::log(spec.upper - spec.lower) - std::log1p(std::exp(-z)) - std::log1p(std::exp(z));
}

double ElementwiseTransformComponent::dx_dz(const Spec &spec, double z)
{
  switch (spec.kind)
  {
    case Kind::log:
    case Kind::lower: return std::exp(z);
    case Kind::upper: return -std::exp(z);
    case Kind::logit:
    {
      double s = 1.0/(1.0 + std::exp(-z));
      return (spec.upper - spec.lower)*s*(1.0 - s);
    }
  }
  return 1.0;
}

double ElementwiseTransformComponent::d_log_dx_dz(const Spec &spec, double z)
{
  if (spec.kind != Kind::logit)
    return 1.0;
  return 1.0 - 2.0/(1.0 + std::exp(-z));
}

Parameters ElementwiseTransformComponent::transform(const Parameters &input) const
{
  // new matrices for the transformed variables (the input's are shared)
  Parameters changed;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* x = input.find_vector(spec.variable);
    if (x == NULL)
      continue;
    arma::mat z(arma::size(*x));
    for (arma::uword k = 0; k < x->n_elem; ++k)
      z[k] = forward(spec, (*x)[k]);
    changed[spec.variable] = z;
  }
  Parameters output = input;
  output.deep_overwrite_with_variables_in_argument(changed);
  return output;
}

Parameters ElementwiseTransformComponent::inverse_transform(const Parameters &input) const
{
  Parameters changed;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* z = input.find_vector(spec.variable);
    if (z == NULL)
      continue;
    arma::mat x(arma::size(*z));
    for (arma::uword k = 0; k < z->n_elem; ++k)
      x[k] = backward(spec, (*z)[k]);
    changed[spec.variable] = x;
  }
  Parameters output = input;
  output.deep_overwrite_with_variables_in_argument(changed);
  return output;
}

arma::mat ElementwiseTransformComponent::jacobian(const Parameters &input) const
{
  // dz/dx on the transformed variables' elements (input: the model's scale)
  std::vector<double> diagonal;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* x = input.find_vector(spec.variable);
    if (x != NULL)
      for (arma::uword k = 0; k < x->n_elem; ++k)
        diagonal.push_back(1.0/dx_dz(spec, forward(spec, (*x)[k])));
  }
  return arma::diagmat(arma::colvec(diagonal));
}

arma::mat ElementwiseTransformComponent::inverse_jacobian(const Parameters &input) const
{
  // dx/dz on the transformed variables' elements (input: the transformed scale)
  std::vector<double> diagonal;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* z = input.find_vector(spec.variable);
    if (z != NULL)
      for (arma::uword k = 0; k < z->n_elem; ++k)
        diagonal.push_back(dx_dz(spec, (*z)[k]));
  }
  return arma::diagmat(arma::colvec(diagonal));
}

double ElementwiseTransformComponent::log_abs_jacobian_determinant(const Parameters &input) const
{
  double output = 0.0;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* x = input.find_vector(spec.variable);
    if (x != NULL)
      for (arma::uword k = 0; k < x->n_elem; ++k)
        output = output - log_dx_dz(spec, forward(spec, (*x)[k]));
  }
  return output;
}

double ElementwiseTransformComponent::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  double output = 0.0;
  for (const Spec &spec : this->specs)
  {
    const arma::mat* z = parameters.find_vector(spec.variable);
    if (z != NULL)
      for (arma::uword k = 0; k < z->n_elem; ++k)
        output = output + log_dx_dz(spec, (*z)[k]);
  }
  return output;
}

void ElementwiseTransformComponent::gradient_to_transformed(const std::string &variable,
                                                            const Parameters &transformed,
                                                            arma::mat &gradient) const
{
  const Spec* spec = this->find(variable);
  if (spec == NULL)
    return; // (identity)
  const arma::mat* z = transformed.find_vector(variable);
  if ((z == NULL) || (z->n_elem != gradient.n_elem))
    Rcpp::stop("The transform of " + variable + " does not match its gradient.");
  for (arma::uword k = 0; k < gradient.n_elem; ++k)
    gradient[k] = gradient[k]*dx_dz(*spec, (*z)[k]) + d_log_dx_dz(*spec, (*z)[k]);
}

int ElementwiseTransformComponent::find_child_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

int ElementwiseTransformComponent::find_parent_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

void ElementwiseTransformComponent::set_child(TransformComponent* new_child)
{
}

void ElementwiseTransformComponent::set_parent(TransformComponent* new_parent)
{
}

}
