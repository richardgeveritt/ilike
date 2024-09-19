// Include the class definition.
#include "custom_root_transform_component.h"

namespace ilike
{
// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
CustomRootTransformComponent::CustomRootTransformComponent()
:TransformComponent()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
CustomRootTransformComponent::~CustomRootTransformComponent()
{
  // Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
CustomRootTransformComponent::CustomRootTransformComponent(const CustomRootTransformComponent &another)
{
  this->make_copy(another);
}

// Returns a pointer to a base class of type CustomRootTransformComponent.
CustomRootTransformComponent* CustomRootTransformComponent::duplicate() const
{
  return new CustomRootTransformComponent(*this);
}

// Copy all the members of the class.
void CustomRootTransformComponent::make_copy(const CustomRootTransformComponent &another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->transform_function = another.transform_function;
  this->jacobian_function = another.jacobian_function;
  this->inverse_transform_function = another.inverse_transform_function;
  this->inverse_jacobian_function = another.inverse_jacobian_function;
  
  // child dealt with in Transform
}

// The = operator.
void CustomRootTransformComponent::operator=(const CustomRootTransformComponent &another)
{
  if(this==&another)//a==a
    return;
  
  TransformComponent::operator=(another);
  
  this->make_copy(another);
}

Parameters CustomRootTransformComponent::transform(const Parameters &parameters) const
{
  return this->transform_function(this->child->transform(parameters));
}

arma::mat CustomRootTransformComponent::jacobian(const Parameters &parameters) const
{
  return this->jacobian_function(this->child->transform(parameters))*this->child->jacobian(parameters);
}

double CustomRootTransformComponent::log_abs_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->jacobian_function(this->child->transform(parameters))).real() + this->child->log_abs_jacobian_determinant(parameters);
}

Parameters CustomRootTransformComponent::inverse_transform(const Parameters &parameters) const
{
  return this->inverse_transform_function(parameters);
}

arma::mat CustomRootTransformComponent::inverse_jacobian(const Parameters &parameters) const
{
  return this->inverse_jacobian_function(parameters);
}

double CustomRootTransformComponent::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->inverse_jacobian_function(parameters)).real();
}

int CustomRootTransformComponent::find_child_index(const std::vector<TransformComponent*> &components) const
{
  for (size_t i=0; i<components.size(); ++i)
  {
    if (this->child==components[i])
      return i;
  }
  return -1;
}

int CustomRootTransformComponent::find_parent_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

void CustomRootTransformComponent::set_child(TransformComponent* new_child)
{
  this->child = new_child;
}

void CustomRootTransformComponent::set_parent(TransformComponent* new_parent)
{
  Rcpp::stop("CustomLeafRootTransformComponent::set_parent - this component has no children.");
}
}
