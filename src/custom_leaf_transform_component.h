#ifndef CUSTOMLEAFTRANSFORMCOMPONENT_H
#define CUSTOMLEAFTRANSFORMCOMPONENT_H

#include <vector>
#include "transform_component.h"
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file custom_leaf_transform_component.h
   * @brief Defines the CustomLeafTransformComponent class.
   *
   * A custom leaf transform component. Applies a custom leaf transformation to parameters as part of a composite Transform.
   *
   * @namespace ilike
   * @class CustomLeafTransformComponent
   * @brief A custom leaf transform component derived from TransformComponent.
   */


class CustomLeafTransformComponent : public TransformComponent
{
public:
  
  /**
   * @brief Default constructor for CustomLeafTransformComponent.
   */
  CustomLeafTransformComponent();
  /**
   * @brief Destructor for CustomLeafTransformComponent.
   */
  virtual ~CustomLeafTransformComponent();
  
  /**
   * @brief Copy constructor for CustomLeafTransformComponent.
   *
   * @param another The CustomLeafTransformComponent instance to copy from.
   */
  CustomLeafTransformComponent(const CustomLeafTransformComponent &another);
  /**
   * @brief Creates a deep copy of this CustomLeafTransformComponent object.
   *
   * @return The result.
   */
  CustomLeafTransformComponent* duplicate() const;
  /**
   * @brief Copies the state of another CustomLeafTransformComponent into this object.
   *
   * @param another The CustomLeafTransformComponent instance to copy from.
   */
  void make_copy(const CustomLeafTransformComponent &another);
  /**
   * @brief Assignment operator for CustomLeafTransformComponent.
   *
   * @param another The CustomLeafTransformComponent instance to copy from.
   */
  void operator=(const CustomLeafTransformComponent &another);
  
  /**
   * @brief Performs the transform operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  Parameters transform(const Parameters &input) const;
  /**
   * @brief Performs the inverse transform operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  Parameters inverse_transform(const Parameters &input) const;
  /**
   * @brief Performs the jacobian operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  arma::mat jacobian(const Parameters &input) const;
  /**
   * @brief Performs the inverse jacobian operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  arma::mat inverse_jacobian(const Parameters &input) const;
  /**
   * @brief Performs the log abs jacobian determinant operation.
   *
   * @param input The input.
   *
   * @return The result.
   */
  double log_abs_jacobian_determinant(const Parameters &input) const;
  /**
   * @brief Performs the log abs inverse jacobian determinant operation.
   *
   * @param parameters The parameters.
   *
   * @return The result.
   */
  double log_abs_inverse_jacobian_determinant(const Parameters &parameters) const;
  
  /**
   * @brief Finds child index.
   *
   * @param components The components.
   *
   * @return The result.
   */
  int find_child_index(const std::vector<TransformComponent*> &components) const;
  /**
   * @brief Finds parent index.
   *
   * @param components The components.
   *
   * @return The result.
   */
  int find_parent_index(const std::vector<TransformComponent*> &components) const;
  /**
   * @brief Sets the child.
   *
   * @param new_child The new child.
   */
  void set_child(TransformComponent* new_child);
  /**
   * @brief Sets the parent.
   *
   * @param new_parent The new parent.
   */
  void set_parent(TransformComponent* new_parent);
  
protected: // Things that can be accessed in this class and subclasses.
  
  /** @brief The transform function. */
  TransformPtr transform_function;
  /** @brief The jacobian function. */
  JacobianPtr jacobian_function;
  
  /** @brief The inverse transform function. */
  TransformPtr inverse_transform_function;
  /** @brief The inverse jacobian function. */
  JacobianPtr inverse_jacobian_function;
  
  /** @brief The parent. */
  TransformComponent* parent;
  
};

}

#endif
