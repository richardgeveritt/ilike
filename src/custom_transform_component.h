#ifndef CUSTOMTRANSFORMCOMPONENT_H
#define CUSTOMTRANSFORMCOMPONENT_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "transform_component.h"
#include "ilike_header.h"

class CustomTransformComponent : public TransformComponent
{
public:
  
  CustomTransformComponent();
  
  virtual ~CustomTransformComponent();
  
  CustomTransformComponent(const CustomTransformComponent &another);
  CustomTransformComponent* duplicate() const;
  void make_copy(const CustomTransformComponent &another);
  void operator=(const CustomTransformComponent &another);
  
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
  
  TransformComponent* parent;
  TransformComponent* child;
  
};

#endif
