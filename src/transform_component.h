#ifndef TRANSFORMCOMPONENT_H
#define TRANSFORMCOMPONENT_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file transform_component.h
   * @brief Defines the TransformComponent class.
   *
   * A generic transform component. Applies a generic transformation to parameters as part of a composite Transform.
   *
   * @namespace ilike
   * @class TransformComponent
   * @brief The transform component class.
   */


class TransformComponent
{
public:
  
  /**
   * @brief Default constructor for TransformComponent.
   */
  TransformComponent();
  /**
   * @brief Destructor for TransformComponent.
   */
  virtual ~TransformComponent();
  
  /**
   * @brief Copy constructor for TransformComponent.
   *
   * @param another The TransformComponent instance to copy from.
   */
  TransformComponent(const TransformComponent &another);
  /**
   * @brief Creates a deep copy of this TransformComponent object.
   *
   * @return The result.
   */
  virtual TransformComponent* duplicate() const=0;
  /**
   * @brief Copies the state of another TransformComponent into this object.
   *
   * @param another The TransformComponent instance to copy from.
   */
  void make_copy(const TransformComponent &another);
  /**
   * @brief Assignment operator for TransformComponent.
   *
   * @param another The TransformComponent instance to copy from.
   */
  void operator=(const TransformComponent &another);
  
  /**
   * @brief Performs the transform operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual Parameters transform(const Parameters &input) const=0;
  /**
   * @brief Performs the inverse transform operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual Parameters inverse_transform(const Parameters &input) const=0;
  /**
   * @brief Performs the jacobian operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual arma::mat jacobian(const Parameters &input) const=0;
  /**
   * @brief Performs the inverse jacobian operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual arma::mat inverse_jacobian(const Parameters &input) const=0;
  /**
   * @brief Performs the log abs jacobian determinant operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  virtual double log_abs_jacobian_determinant(const Parameters &input) const=0;
  /**
   * @brief Performs the log abs inverse jacobian determinant operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  virtual double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const=0;
  
  /**
   * @brief Finds child index.
   *
   * @param components The components.
   *
   * @return The result.
   */
  virtual int find_child_index(const std::vector<TransformComponent*> &components) const=0;
  /**
   * @brief Finds parent index.
   *
   * @param components The components.
   *
   * @return The result.
   */
  virtual int find_parent_index(const std::vector<TransformComponent*> &components) const=0;
  /**
   * @brief Sets the child.
   *
   * @param new_child The new child.
   */
  virtual void set_child(TransformComponent* new_child)=0;
  /**
   * @brief Sets the parent.
   *
   * @param new_parent The new parent.
   */
  virtual void set_parent(TransformComponent* new_parent)=0;
  
protected: // Things that can be accessed in this class and subclasses.
  
};
}

#endif
