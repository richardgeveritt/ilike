#ifndef TRANSFORMCOMPONENT_H
#define TRANSFORMCOMPONENT_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "ilike_header.h"

namespace ilike
{
class TransformComponent
{
public:
  
  TransformComponent();
  virtual ~TransformComponent();
  
  TransformComponent(const TransformComponent &another);
  virtual TransformComponent* duplicate() const=0;
  void make_copy(const TransformComponent &another);
  void operator=(const TransformComponent &another);
  
  virtual Parameters transform(const Parameters &input) const=0;
  virtual Parameters inverse_transform(const Parameters &input) const=0;
  virtual arma::mat jacobian(const Parameters &input) const=0;
  virtual arma::mat inverse_jacobian(const Parameters &input) const=0;
  virtual double log_abs_jacobian_determinant(const Parameters &input) const=0;
  virtual double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const=0;
  
  virtual int find_child_index(const std::vector<TransformComponent*> &components) const=0;
  virtual int find_parent_index(const std::vector<TransformComponent*> &components) const=0;
  virtual void set_child(TransformComponent* new_child)=0;
  virtual void set_parent(TransformComponent* new_parent)=0;
  
protected: // Things that can be accessed in this class and subclasses.
  
};
}

#endif
