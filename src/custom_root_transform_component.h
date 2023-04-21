#ifndef CUSTOMROOTTRANSFORMCOMPONENT_H
#define CUSTOMROOTTRANSFORMCOMPONENT_H

#include <vector>
#include "transform_component.h"
#include "ilike_header.h"

class CustomRootTransformComponent : public TransformComponent
{
public:
  
  CustomRootTransformComponent();
  virtual ~CustomRootTransformComponent();
  
  CustomRootTransformComponent(const CustomRootTransformComponent &another);
  CustomRootTransformComponent* duplicate() const;
  void make_copy(const CustomRootTransformComponent &another);
  void operator=(const CustomRootTransformComponent &another);
  
  Parameters transform(const Parameters &input) const;
  Parameters inverse_transform(const Parameters &input) const;
  arma::mat jacobian(const Parameters &input) const;
  arma::mat inverse_jacobian(const Parameters &input) const;
  double log_abs_jacobian_determinant(const Parameters &input) const;
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;
  
  int find_child_index(const std::vector<TransformComponent*> &components) const;
  int find_parent_index(const std::vector<TransformComponent*> &components) const;
  void set_child(TransformComponent* new_child);
  void set_parent(TransformComponent* new_parent);
  
protected: // Things that can be accessed in this class and subclasses.
  
  TransformPtr transform_function;
  JacobianPtr jacobian_function;
  
  TransformPtr inverse_transform_function;
  JacobianPtr inverse_jacobian_function;
  
  TransformComponent* child;
  
};


#endif
