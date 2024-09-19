// Include the class definition.
#include "custom_leaf_transform_component.h"

namespace ilike
{
// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
CustomLeafTransformComponent::CustomLeafTransformComponent()
:TransformComponent()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
CustomLeafTransformComponent::~CustomLeafTransformComponent()
{
  // Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
CustomLeafTransformComponent::CustomLeafTransformComponent(const CustomLeafTransformComponent &another)
{
  this->make_copy(another);
}

// Returns a pointer to a base class of type CustomLeafTransformComponent.
CustomLeafTransformComponent* CustomLeafTransformComponent::duplicate() const
{
  return new CustomLeafTransformComponent(*this);
}

// Copy all the members of the class.
void CustomLeafTransformComponent::make_copy(const CustomLeafTransformComponent &another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->transform_function = another.transform_function;
  this->jacobian_function = another.jacobian_function;
  this->inverse_transform_function = another.inverse_transform_function;
  this->inverse_jacobian_function = another.inverse_jacobian_function;
  
  // parent dealt with in Transform
}

// The = operator.
void CustomLeafTransformComponent::operator=(const CustomLeafTransformComponent &another)
{
  if(this==&another)//a==a
    return;
  
  TransformComponent::operator=(another);
  
  this->make_copy(another);
}

Parameters CustomLeafTransformComponent::transform(const Parameters &parameters) const
{
  return this->transform_function(parameters);
}

arma::mat CustomLeafTransformComponent::jacobian(const Parameters &parameters) const
{
  return this->jacobian_function(parameters);
}

double CustomLeafTransformComponent::log_abs_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->jacobian_function(parameters)).real();
}

Parameters CustomLeafTransformComponent::inverse_transform(const Parameters &parameters) const
{
  return this->inverse_transform_function(this->parent->inverse_transform(parameters));
}

arma::mat CustomLeafTransformComponent::inverse_jacobian(const Parameters &parameters) const
{
  return this->inverse_jacobian_function(this->parent->inverse_transform(parameters))*this->parent->inverse_jacobian(parameters);
}

double CustomLeafTransformComponent::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->inverse_jacobian_function(this->parent->inverse_transform(parameters))).real() + this->parent->log_abs_inverse_jacobian_determinant(parameters);
}

int CustomLeafTransformComponent::find_child_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

int CustomLeafTransformComponent::find_parent_index(const std::vector<TransformComponent*> &components) const
{
  for (size_t i=0; i<components.size(); ++i)
  {
    if (this->parent==components[i])
      return i;
  }
  return -1;
}

void CustomLeafTransformComponent::set_child(TransformComponent* new_child)
{
  Rcpp::stop("CustomLeafRootTransformComponent::set_child - this component has no children.");
}

void CustomLeafTransformComponent::set_parent(TransformComponent* new_parent)
{
  this->parent = new_parent;
}
}
