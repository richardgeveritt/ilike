#ifndef ELEMENTWISETRANSFORMCOMPONENT_H
#define ELEMENTWISETRANSFORMCOMPONENT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>

#include "transform_component.h"

namespace ilike
{

// The built-in constraining transforms (dev/plan.md 3.18 H2, C209), applied
// element by element to named variables (the others are unchanged):
// - log:          z = log(x)                    x > 0
// - lower(a):     z = log(x - a)                x > a
// - upper(b):     z = log(b - x)                x < b
// - logit(a, b):  z = logit((x - a)/(b - a))    a < x < b
// The Jacobian is diagonal, so a gradient on the model's scale gives the
// gradient on the transformed scale exactly (gradient_to_transformed).
// And three applied to a whole variable (C232), as Stan's:
// - simplex:  K probabilities summing to 1 <- K - 1 values by stick-breaking,
//             y_k = logit^-1(z_k - log(K - k)), x_k = y_k (1 - sum_{j<k} x_j)
// - cov:      a p x p positive definite matrix <- p(p + 1)/2 values: its
//             Cholesky factor L, log on the diagonal (column by column)
// - corr:     a p x p correlation matrix <- p(p - 1)/2 values: canonical
//             partial correlations tanh(z) (column by column), turned into a
//             Cholesky factor (Stan's read_corr_L)
// Their log |dx/dz| are Stan's (with respect to the free elements: the first
// K - 1 probabilities; the lower triangle; the strict lower triangle). Their
// gradient's chain rule uses dx/dz and the gradient of log |dx/dz| by
// central differences (which affects only the efficiency of gradient-based
// moves: their exactness rests on the density, not on the gradient).
// jacobian()/inverse_jacobian() are the diagonal matrices over the transformed
// variables' elements, in the order of the specification.
class ElementwiseTransformComponent : public TransformComponent
{
public:

  enum class Kind { log, lower, upper, logit, simplex, cov, corr };

  struct Spec
  {
    std::string variable;
    Kind kind;
    double lower;
    double upper;
  };

  ElementwiseTransformComponent();
  ElementwiseTransformComponent(const std::vector<Spec> &specs_in);
  virtual ~ElementwiseTransformComponent();

  ElementwiseTransformComponent(const ElementwiseTransformComponent &another);
  TransformComponent* duplicate() const;

  Parameters transform(const Parameters &input) const;
  Parameters inverse_transform(const Parameters &input) const;
  arma::mat jacobian(const Parameters &input) const;
  arma::mat inverse_jacobian(const Parameters &input) const;
  // log |dz/dx| at x (input on the model's scale)
  double log_abs_jacobian_determinant(const Parameters &input) const;
  // log |dx/dz| at z (parameters on the transformed scale)
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;

  int find_child_index(const std::vector<TransformComponent*> &components) const;
  int find_parent_index(const std::vector<TransformComponent*> &components) const;
  void set_child(TransformComponent* new_child);
  void set_parent(TransformComponent* new_parent);

  // the gradient of log pi(x(z)) + log |dx/dz| with respect to `variable` on
  // the transformed scale, from the gradient of log pi with respect to x
  // (`gradient`, replaced), at the transformed point `transformed`
  void gradient_to_transformed(const std::string &variable,
                               const Parameters &transformed,
                               arma::mat &gradient) const;

  const std::vector<Spec> &get_specs() const { return this->specs; }
  // the number of elements of `variable` on the model's scale, from its
  // number on the transformed scale (they differ for simplex, cov, corr)
  arma::uword model_size(const std::string &variable, arma::uword transformed_size) const;

private:

  const Spec* find(const std::string &variable) const;
  static double forward(const Spec &spec, double x);
  static double backward(const Spec &spec, double z);
  static double log_dx_dz(const Spec &spec, double z);       // log |dx/dz|
  static double dx_dz(const Spec &spec, double z);
  static double d_log_dx_dz(const Spec &spec, double z);     // d/dz log |dx/dz|
  static bool whole(Kind kind) { return (kind == Kind::simplex) || (kind == Kind::cov) || (kind == Kind::corr); }
  static arma::colvec whole_forward(const Spec &spec, const arma::mat &x);
  static arma::mat whole_backward(const Spec &spec, const arma::colvec &z);
  static double whole_log_dx_dz(const Spec &spec, const arma::colvec &z);

  std::vector<Spec> specs;
};

}

#endif
