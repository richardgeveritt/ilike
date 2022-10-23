// Include the class definition.
#include "transform.h"

// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Transform::Transform(void)
{
	// Set all the pointers owned by this class to be something (possibly NULL).
}

// Comment about function.
Transform::~Transform(void)
{
	// Delete all the pointers owned by this class.
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
Transform::Transform(const Transform &another)
{
	make_copy(another);
}

// Returns a pointer to a base class of type Transform.
Transform* Transform::duplicate() const
{
	return new Transform(*this);
}

// Copy all the members of the class.
void Transform::make_copy(const Transform &another)
{
	// Copy all members, duplicating the memory where appropriate.
  this->transform = another.transform;
  this->jacobian = another.jacobian;
  this->inverse_transform = another.inverse_transform;
  this->inverse_jacobian = another.inverse_jacobian;
}

// The = operator.
void Transform::operator=(const Transform &another)
{
	if(this==&another)//a==a
		return;
	
	make_copy(another);
}

/*
Parameters Transform::transform(const Parameters &parameters) const
{
  
}

arma::mat Transform::jacobian(const Parameters &parameters) const
{
  
}
*/

double Transform::log_abs_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->jacobian(parameters)).real();
}

/*
Parameters Transform::inverse_transform(const Parameters &parameters) const
{
  
}

arma::mat inverse_jacobian(const Parameters &parameters) const
{
  
}
*/

double Transform::log_abs_inverse_jacobian_determinant(const Parameters &parameters) const
{
  return log_det(this->inverse_jacobian(parameters)).real();
}
