// Include the class definition.
#include "transform.h"
#include "transform_component.h"
#include "custom_transform_component.h"
#include "custom_leafroot_transform_component.h"

namespace ilike
{
// Include other class definitions that are needed.

// Public functions //

// Default constructor and destructor.

// Comment about function.
Transform::Transform()
{
  // Set all the pointers owned by this class to be something (possibly NULL).
}

Transform::Transform(TransformPtr transform_in)
{
  this->components.emplace_back(new CustomLeafRootTransformComponent(transform_in));
  this->root = this->components.begin()->get();
  this->leaf = this->components.begin()->get();
}

Transform::Transform(TransformPtr transform_in,
                     TransformPtr inverse_transform_in)
{
  this->components.emplace_back(new CustomLeafRootTransformComponent(transform_in,
                                                                  inverse_transform_in));
  this->root = this->components.begin()->get();
  this->leaf = this->components.begin()->get();
}

// Comment about function.
Transform::~Transform()
{
  // Delete all the pointers owned by this class.
  
  for (auto i=this->components.begin();
       i!=this->components.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  
}

// Other constuctors.

// Everything you need to copy the class.

// The copy constructor.
Transform::Transform(const Transform &another)
{
  this->make_copy(another);
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
  
  this->components.resize(0);
  this->components.reserve(another.components.size());
  
  std::vector<int> child_indices;
  child_indices.reserve(another.components.size());
  
  std::vector<int> parent_indices;
  parent_indices.reserve(another.components.size());

  // the components' own lookup functions take plain (non-owning) pointers
  std::vector<TransformComponent*> another_components;
  another_components.reserve(another.components.size());
  for (const auto &component : another.components)
    another_components.push_back(component.get());
  
  for (size_t i=0;
       i<another.components.size();
       ++i)
  {
    if (another.components[i]!=NULL)
      this->components.emplace_back(another.components[i]->duplicate());
    else
      this->components.emplace_back();
    
    child_indices.push_back(another.components[i]->find_child_index(another_components));
    parent_indices.push_back(another.components[i]->find_parent_index(another_components));
    
    if (another.components[i].get()==another.root)
      this->root = this->components.back().get();
    
    if (another.components[i].get()==another.leaf)
      this->leaf = this->components.back().get();
  }
  
  for (size_t i=0;
       i<another.components.size();
       ++i)
  {
    if (child_indices[i]>=0)
      this->components[i]->set_child(this->components[child_indices[i]].get());
    
    if (parent_indices[i]>=0)
      this->components[i]->set_parent(this->components[parent_indices[i]].get());
  }
}

// The = operator.
void Transform::operator=(const Transform &another)
{
  if(this==&another)//a==a
    return;
  
  for (auto i=this->components.begin();
       i!=this->components.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  this->components.clear();
  
  this->make_copy(another);
}

Parameters Transform::transform(const Parameters &input) const
{
  return this->root->transform(input);
}

Parameters Transform::inverse_transform(const Parameters &input) const
{
  return this->leaf->inverse_transform(input);
}

arma::mat Transform::jacobian(const Parameters &input) const
{
  return this->root->jacobian(input);
}

arma::mat Transform::inverse_jacobian(const Parameters &input) const
{
  return this->leaf->jacobian(input);
}

double Transform::log_abs_jacobian_determinant(const Parameters &input) const
{
  return this->root->log_abs_jacobian_determinant(input);
}

double Transform::log_abs_inverse_jacobian_determinant(const Parameters &input) const
{
  return this->leaf->log_abs_inverse_jacobian_determinant(input);
}
}
