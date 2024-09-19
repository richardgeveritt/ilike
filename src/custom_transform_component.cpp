// Include the class definition.
#include "custom_transform_component.h"

namespace ilike
{
// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
CustomTransformComponent::CustomTransformComponent()
:TransformComponent()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
CustomTransformComponent::~CustomTransformComponent()
{
  // Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
CustomTransformComponent::CustomTransformComponent(const CustomTransformComponent &another)
{
  this->make_copy(another);
}

// Returns a pointer to a base class of type CustomTransformComponent.
CustomTransformComponent* CustomTransformComponent::duplicate() const
{
  return new CustomTransformComponent(*this);
}

// Copy all the members of the class.
void CustomTransformComponent::make_copy(const CustomTransformComponent &another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->transform_function = another.transform_function;
  this->jacobian_function = another.jacobian_function;
  this->inverse_transform_function = another.inverse_transform_function;
  this->inverse_jacobian_function = another.inverse_jacobian_function;
  
  // parent and child dealt with in Transform
}

// The = operator.
void CustomTransformComponent::operator=(const CustomTransformComponent &another)
{
  if(this==&another)//a==a
    return;
  
  TransformComponent::operator=(another);
  
  this->make_copy(another);
}

Parameters CustomTransformComponent::transform(const Parameters &parameters) const
{
  return this->transform_function(this->child->transform(parameters));
}

arma::mat CustomTransformComponent::jacobian(const Parameters &parameters) const
{
  return this->jacobian_function(this->child->transform(parameters))*this->child->jacobian(parameters);
}

double CustomTransformComponent::log_abs_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->jacobian_function(this->child->transform(parameters))).real() + this->child->log_abs_jacobian_determinant(parameters);
}

Parameters CustomTransformComponent::inverse_transform(const Parameters &parameters) const
{
  return this->inverse_transform_function(this->parent->inverse_transform(parameters));
}

arma::mat CustomTransformComponent::inverse_jacobian(const Parameters &parameters) const
{
  return this->inverse_jacobian_function(this->parent->inverse_transform(parameters))*this->parent->inverse_jacobian(parameters);
}

double CustomTransformComponent::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->inverse_jacobian_function(this->parent->inverse_transform(parameters))).real() + this->parent->log_abs_inverse_jacobian_determinant(parameters);
}

int CustomTransformComponent::find_child_index(const std::vector<TransformComponent*> &components) const
{
  for (size_t i=0; i<components.size(); ++i)
  {
    if (this->child==components[i])
      return i;
  }
  return -1;
}

int CustomTransformComponent::find_parent_index(const std::vector<TransformComponent*> &components) const
{
  for (size_t i=0; i<components.size(); ++i)
  {
    if (this->parent==components[i])
      return i;
  }
  return -1;
}

void CustomTransformComponent::set_child(TransformComponent* new_child)
{
  this->child = new_child;
}

void CustomTransformComponent::set_parent(TransformComponent* new_parent)
{
  this->parent = new_parent;
}
}
