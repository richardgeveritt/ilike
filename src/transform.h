#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "ilike_header.h"

class TransformComponent;

class Transform
{
public:

	Transform();
  
  Transform(TransformPtr transform_in);
  Transform(TransformPtr transform_in,
            TransformPtr inverse_transform_in);
	virtual ~Transform();

	Transform(const Transform &another);
	Transform* duplicate() const;
	void make_copy(const Transform &another);
	void operator=(const Transform &another);
  
  Parameters transform(const Parameters &input) const;
  Parameters inverse_transform(const Parameters &input) const;
  arma::mat jacobian(const Parameters &input) const;
  arma::mat inverse_jacobian(const Parameters &input) const;
  double log_abs_jacobian_determinant(const Parameters &input) const;
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;

protected: // Things that can be accessed in this class and subclasses.

  // stored here
  std::vector<TransformComponent*> components;
  
  // points to one of the components
  TransformComponent* root;
  
  // points to one of the components
  TransformComponent* leaf;
  
};

#endif
