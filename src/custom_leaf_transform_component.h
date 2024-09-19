#ifndef CUSTOMLEAFTRANSFORMCOMPONENT_H
#define CUSTOMLEAFTRANSFORMCOMPONENT_H

#include <vector>
#include "transform_component.h"
#include "ilike_header.h"

namespace ilike
{
class CustomLeafTransformComponent : public TransformComponent
{
public:
  
  CustomLeafTransformComponent();
  virtual ~CustomLeafTransformComponent();
  
  CustomLeafTransformComponent(const CustomLeafTransformComponent &another);
  CustomLeafTransformComponent* duplicate() const;
  void make_copy(const CustomLeafTransformComponent &another);
  void operator=(const CustomLeafTransformComponent &another);
  
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
  
};

}

#endif
