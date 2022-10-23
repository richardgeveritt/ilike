#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "function_pointers.h"

class Transform
{
public:

	Transform();
	virtual ~Transform();

	Transform(const Transform &another);
	Transform* duplicate() const;
	void make_copy(const Transform &another);
	void operator=(const Transform &another);
  
  //Parameters transform(const Parameters &parameters) const;
  //arma::mat jacobian(const Parameters &parameters) const;
  TransformPtr transform;
  JacobianPtr jacobian;
  double log_abs_jacobian_determinant(const Parameters &parameters) const;
  
  //Parameters inverse_transform(const Parameters &parameters) const;
  //arma::mat inverse_jacobian(const Parameters &parameters) const;
  TransformPtr inverse_transform;
  JacobianPtr inverse_jacobian;
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;

protected: // Things that can be accessed in this class and subclasses.

};

#endif
