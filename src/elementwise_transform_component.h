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
// jacobian()/inverse_jacobian() are the diagonal matrices over the transformed
// variables' elements, in the order of the specification.
class ElementwiseTransformComponent : public TransformComponent
{
public:

  enum class Kind { log, lower, upper, logit };

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

private:

  const Spec* find(const std::string &variable) const;
  static double forward(const Spec &spec, double x);
  static double backward(const Spec &spec, double z);
  static double log_dx_dz(const Spec &spec, double z);       // log |dx/dz|
  static double dx_dz(const Spec &spec, double z);
  static double d_log_dx_dz(const Spec &spec, double z);     // d/dz log |dx/dz|

  std::vector<Spec> specs;
};

}

#endif
