#ifndef CUSTOMLEAFROOTTRANSFORMCOMPONENT_H
#define CUSTOMLEAFROOTTRANSFORMCOMPONENT_H

#include <vector>
#include "transform_component.h"
#include "ilike_header.h"

namespace ilike
{
class CustomLeafRootTransformComponent : public TransformComponent
{
public:
  
  CustomLeafRootTransformComponent();
  
  CustomLeafRootTransformComponent(TransformPtr transform_function_in);
  CustomLeafRootTransformComponent(TransformPtr transform_function_in,
                                   TransformPtr inverse_transform_function_in);
  
  virtual ~CustomLeafRootTransformComponent();
  
  CustomLeafRootTransformComponent(const CustomLeafRootTransformComponent &another);
  CustomLeafRootTransformComponent* duplicate() const;
  void make_copy(const CustomLeafRootTransformComponent &another);
  void operator=(const CustomLeafRootTransformComponent &another);
  
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
  
};
}


#endif
