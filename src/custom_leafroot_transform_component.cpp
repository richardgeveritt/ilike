// Include the class definition.
#include "custom_leafroot_transform_component.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
CustomLeafRootTransformComponent::CustomLeafRootTransformComponent()
:TransformComponent()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

CustomLeafRootTransformComponent::CustomLeafRootTransformComponent(TransformPtr transform_function_in)
:TransformComponent()
{
  this->transform_function = transform_function_in;
}

CustomLeafRootTransformComponent::CustomLeafRootTransformComponent(TransformPtr transform_function_in,
                                                                   TransformPtr inverse_transform_function_in)
:TransformComponent()
{
  this->transform_function = transform_function_in;
  this->inverse_transform_function = inverse_transform_function_in;
}

// Comment about function.
CustomLeafRootTransformComponent::~CustomLeafRootTransformComponent()
{
  // Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
CustomLeafRootTransformComponent::CustomLeafRootTransformComponent(const CustomLeafRootTransformComponent &another)
{
  this->make_copy(another);
}

// Returns a pointer to a base class of type CustomLeafRootTransformComponent.
CustomLeafRootTransformComponent* CustomLeafRootTransformComponent::duplicate() const
{
  return new CustomLeafRootTransformComponent(*this);
}

// Copy all the members of the class.
void CustomLeafRootTransformComponent::make_copy(const CustomLeafRootTransformComponent &another)
{
  // Copy all members, duplicating the memory where appropriate.
  this->transform_function = another.transform_function;
  this->jacobian_function = another.jacobian_function;
  this->inverse_transform_function = another.inverse_transform_function;
  this->inverse_jacobian_function = another.inverse_jacobian_function;
  
  // child dealt with in Transform
}

// The = operator.
void CustomLeafRootTransformComponent::operator=(const CustomLeafRootTransformComponent &another)
{
  if(this==&another)//a==a
    return;
  
  TransformComponent::operator=(another);
  
  this->make_copy(another);
}

Parameters CustomLeafRootTransformComponent::transform(const Parameters &parameters) const
{
  return this->transform_function(parameters);
}

arma::mat CustomLeafRootTransformComponent::jacobian(const Parameters &parameters) const
{
  return this->jacobian_function(parameters);
}

double CustomLeafRootTransformComponent::log_abs_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->jacobian_function(parameters)).real();
}

Parameters CustomLeafRootTransformComponent::inverse_transform(const Parameters &parameters) const
{
  return this->inverse_transform_function(parameters);
}

arma::mat CustomLeafRootTransformComponent::inverse_jacobian(const Parameters &parameters) const
{
  return this->inverse_jacobian_function(parameters);
}

double CustomLeafRootTransformComponent::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->inverse_jacobian_function(parameters)).real();
}

int CustomLeafRootTransformComponent::find_child_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

int CustomLeafRootTransformComponent::find_parent_index(const std::vector<TransformComponent*> &components) const
{
  return -1;
}

void CustomLeafRootTransformComponent::set_child(TransformComponent* new_child)
{
  Rcpp::stop("CustomLeafRootTransformComponent::set_child - this component has no children.");
}

void CustomLeafRootTransformComponent::set_parent(TransformComponent* new_parent)
{
  Rcpp::stop("CustomLeafRootTransformComponent::set_parent - this component has no children.");
}
