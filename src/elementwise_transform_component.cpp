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

// --- the whole-variable transforms (C232) -------------------------------------

static arma::uword triangle_side(arma::uword n, bool diagonal)
{
  // p with p(p + 1)/2 = n (diagonal) or p(p - 1)/2 = n
  double p = diagonal ? (std::sqrt(8.0*double(n) + 1.0) - 1.0)/2.0 : (std::sqrt(8.0*double(n) + 1.0) + 1.0)/2.0;
  arma::uword side = arma::uword(std::lround(p));
  if ((diagonal ? side*(side + 1)/2 : side*(side - 1)/2) != n)
    Rcpp::stop("A transformed matrix has a number of elements that is not triangular.");
  return side;
}

arma::uword ElementwiseTransformComponent::model_size(const std::string &variable, arma::uword n) const
{
  const Spec* spec = this->find(variable);
  if ((spec == NULL) || !whole(spec->kind))
    return n;
  if (spec->kind == Kind::simplex)
    return n + 1;
  arma::uword p = triangle_side(n, spec->kind == Kind::cov);
  return p*p;
}

arma::colvec ElementwiseTransformComponent::whole_forward(const Spec &spec, const arma::mat &x)
{
  if (spec.kind == Kind::simplex)
  {
    arma::uword K = x.n_elem;
    arma::colvec z(K - 1);
    double stick = 1.0;
    for (arma::uword k = 0; k + 1 < K; ++k)
    {
      double y = x[k]/stick;
      z[k] = std::log(y) - std::log1p(-y) + std::log(double(K - 1 - k));
      stick = stick - x[k];
    }
    return z;
  }
  if (x.n_rows != x.n_cols)
    Rcpp::stop("transform " + std::string(spec.kind == Kind::cov ? "cov" : "corr") + "(" + spec.variable +
               "): the value is not a square matrix.");
  arma::uword p = x.n_rows;
  arma::mat L;
  if (!arma::chol(L, arma::symmatu(x), "lower"))
  {
    // (a matrix that is not positive definite, e.g. at the edge of the support
    // after an extreme step: density 0 there, so no finite transformed value)
    arma::uword n = (spec.kind == Kind::cov) ? p*(p + 1)/2 : p*(p - 1)/2;
    return arma::colvec(n, arma::fill::value(arma::datum::nan));
  }
  std::vector<double> z;
  if (spec.kind == Kind::cov)
  {
    for (arma::uword j = 0; j < p; ++j)
      for (arma::uword i = j; i < p; ++i)
        z.push_back((i == j) ? std::log(L(i, i)) : L(i, j));
    return arma::colvec(z);
  }
  // corr: the canonical partial correlations, column by column
  for (arma::uword j = 0; j + 1 < p; ++j)
    for (arma::uword i = j + 1; i < p; ++i)
    {
      double used = 0.0;
      for (arma::uword k = 0; k < j; ++k)
        used = used + L(i, k)*L(i, k);
      double cpc = L(i, j)/std::sqrt(1.0 - used);
      z.push_back(std::atanh(std::max(-1.0 + 1e-15, std::min(1.0 - 1e-15, cpc))));
    }
  return arma::colvec(z);
}

arma::mat ElementwiseTransformComponent::whole_backward(const Spec &spec, const arma::colvec &z)
{
  if (spec.kind == Kind::simplex)
  {
    arma::uword K = z.n_elem + 1;
    arma::colvec x(K);
    double stick = 1.0;
    for (arma::uword k = 0; k + 1 < K; ++k)
    {
      double y = 1.0/(1.0 + std::exp(-(z[k] - std::log(double(K - 1 - k)))));
      x[k] = stick*y;
      stick = stick - x[k];
    }
    x[K - 1] = stick;
    return x;
  }
  bool covariance = (spec.kind == Kind::cov);
  arma::uword p = triangle_side(z.n_elem, covariance);
  arma::mat L(p, p, arma::fill::zeros);
  if (covariance)
  {
    arma::uword position = 0;
    for (arma::uword j = 0; j < p; ++j)
      for (arma::uword i = j; i < p; ++i)
        L(i, j) = (i == j) ? std::exp(z[position++]) : z[position++];
  }
  else
  {
    // Stan's read_corr_L: column j's entries are cpc * sqrt(what is left of row i's length)
    arma::colvec left(p, arma::fill::ones);
    arma::uword position = 0;
    for (arma::uword j = 0; j < p; ++j)
    {
      L(j, j) = std::sqrt(left[j]);
      for (arma::uword i = j + 1; i < p; ++i)
      {
        double cpc = std::tanh(z[position++]);
        L(i, j) = cpc*std::sqrt(left[i]);
        left[i] = left[i]*(1.0 - cpc*cpc);
      }
    }
  }
  return L*L.t();
}

double ElementwiseTransformComponent::whole_log_dx_dz(const Spec &spec, const arma::colvec &z)
{
  if (spec.kind == Kind::simplex)
  {
    arma::uword K = z.n_elem + 1;
    double stick = 1.0, value = 0.0;
    for (arma::uword k = 0; k + 1 < K; ++k)
    {
      double a = z[k] - std::log(double(K - 1 - k));
      // log y + log(1 - y) + log(stick), y = logit^-1(a)
      value = value - std::log1p(std::exp(-a)) - std::log1p(std::exp(a)) + std::log(stick);
      stick = stick*(1.0 - 1.0/(1.0 + std::exp(-a)));
    }
    return value;
  }
  if (spec.kind == Kind::cov)
  {
    arma::uword p = triangle_side(z.n_elem, true);
    double value = double(p)*std::log(2.0);
    arma::uword position = 0;
    for (arma::uword j = 0; j < p; ++j)
      for (arma::uword i = j; i < p; ++i, ++position)
        if (i == j)
          value = value + double(p - j + 1)*z[position]; // (p - k + 2) log L_kk, k = j + 1, log L_jj = z
    return value;
  }
  // corr: tanh's log(1 - cpc^2), and 0.5 (p - k - 1) log(1 - cpc^2) for column k (1-based)
  arma::uword p = triangle_side(z.n_elem, false);
  double value = 0.0;
  arma::uword position = 0;
  for (arma::uword j = 0; j + 1 < p; ++j)
    for (arma::uword i = j + 1; i < p; ++i, ++position)
    {
      double t = std::tanh(z[position]);
      double log_one_minus = std::log1p(-t*t);
      value = value + log_one_minus + 0.5*double(p - j - 2)*log_one_minus;
    }
  return value;
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
    default: break; // (the whole-variable kinds: whole_forward)
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
    default: break; // (the whole-variable kinds: whole_backward)
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
    default: break; // (the whole-variable kinds: by central differences)
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
    if (whole(spec.kind))
    {
      changed[spec.variable] = arma::mat(whole_forward(spec, *x));
      continue;
    }
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
    if (whole(spec.kind))
    {
      changed[spec.variable] = whole_backward(spec, arma::vectorise(*z));
      continue;
    }
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
    if (whole(spec.kind))
      Rcpp::stop("The Jacobian matrix of the simplex, cov and corr transforms is not available here.");
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
    if (whole(spec.kind))
      Rcpp::stop("The Jacobian matrix of the simplex, cov and corr transforms is not available here.");
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
    if ((x != NULL) && whole(spec.kind))
      output = output - whole_log_dx_dz(spec, whole_forward(spec, *x));
    else if (x != NULL)
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
    if ((z != NULL) && whole(spec.kind))
      output = output + whole_log_dx_dz(spec, arma::vectorise(*z));
    else if (z != NULL)
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
  if (z == NULL)
    Rcpp::stop("The transform of " + variable + " does not match its gradient.");
  if (whole(spec->kind))
  {
    // dx/dz and the gradient of log |dx/dz| by central differences: the
    // gradient on the transformed scale (its size, not the model's)
    arma::colvec z0 = arma::vectorise(*z);
    arma::colvec g = arma::vectorise(gradient);
    if (g.n_elem != this->model_size(variable, z0.n_elem))
      Rcpp::stop("The transform of " + variable + " does not match its gradient.");
    arma::colvec output(z0.n_elem);
    for (arma::uword k = 0; k < z0.n_elem; ++k)
    {
      double h = 1e-6*std::max(1.0, std::abs(z0[k]));
      arma::colvec zp = z0, zm = z0;
      zp[k] = zp[k] + h;
      zm[k] = zm[k] - h;
      arma::colvec dx = (arma::vectorise(whole_backward(*spec, zp)) - arma::vectorise(whole_backward(*spec, zm)))/(2.0*h);
      output[k] = arma::dot(dx, g) + (whole_log_dx_dz(*spec, zp) - whole_log_dx_dz(*spec, zm))/(2.0*h);
    }
    gradient = output;
    return;
  }
  if (z->n_elem != gradient.n_elem)
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
