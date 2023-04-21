// Include the class definition.
#include "transform_component.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
TransformComponent::TransformComponent()
{
	// Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
TransformComponent::~TransformComponent()
{
	// Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
TransformComponent::TransformComponent(const TransformComponent &another)
{
	this->make_copy(another);
}

// Copy all the members of the class.
void TransformComponent::make_copy(const TransformComponent &another)
{
	// Copy all members, duplicating the memory where appropriate.
}

// The = operator.
void TransformComponent::operator=(const TransformComponent &another)
{
	if(this==&another)//a==a
		return;
	
	this->make_copy(another);
}
