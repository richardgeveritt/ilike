#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "ilike_header.h"

namespace ilike
{

class TransformComponent;

class Transform
{
public:
  
  Transform();
  
  Transform(TransformPtr transform_in);
  Transform(TransformPtr transform_in,
            TransformPtr inverse_transform_in);
  // from one component (takes ownership): the built-in elementwise transforms (C209)
  Transform(TransformComponent* component_in);
  virtual ~Transform();
  
  Transform(const Transform &another);
  Transform* duplicate() const;
  void make_copy(const Transform &another);
  void operator=(const Transform &another);
  
  Parameters transform(const Parameters &input) const;
  Parameters inverse_transform(const Parameters &input) const;
  arma::mat jacobian(const Parameters &input) const;
  arma::mat inverse_jacobian(const Parameters &input) const;
  double log_abs_jacobian_determinant(const Parameters &input) const;
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;

  // For the built-in elementwise transforms (C209): the gradient of
  // log pi(x(z)) + log|dx/dz| with respect to `variable` (or a stacked name)
  // on the transformed scale, from the gradient on the model's scale
  // (replaced), at the transformed point. FALSE for other transforms.
  bool gradient_to_transformed(const std::string &variable,
                               const Parameters &transformed,
                               arma::mat &gradient) const;
  // whether this is a built-in elementwise transform
  bool is_elementwise() const;
  
protected: // Things that can be accessed in this class and subclasses.
  
  // stored here
  std::vector<std::unique_ptr<TransformComponent>> components;
  
  // points to one of the components
  TransformComponent* root;
  
  // points to one of the components
  TransformComponent* leaf;
  
};
}

#endif
