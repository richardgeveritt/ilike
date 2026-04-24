#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <vector>
#include <boost/unordered_map.hpp>
#include "ilike_header.h"

namespace ilike
{
  /**
   * @file transform.h
   * @brief Defines the TransformComponent class.
   *
   * A generic transform component. Applies a generic transformation to parameters as part of a composite Transform.
   *
   * @namespace ilike
   * @class TransformComponent
   * @brief The transform component class.
   */


class TransformComponent;

class Transform
{
public:
  
  /**
   * @brief Performs the transform operation.
   */
  Transform();
  
  /**
   * @brief Performs the transform operation.
   *
   * @param transform_in The transform.
   */
  Transform(TransformPtr transform_in);
  Transform(TransformPtr transform_in,
            TransformPtr inverse_transform_in);
  /**
   * @brief Performs the ~transform operation.
   */
  virtual ~Transform();
  
  /**
   * @brief Performs the transform operation.
   *
   * @param another The TransformComponent instance to copy from.
   */
  Transform(const Transform &another);
  /**
   * @brief Creates a deep copy of this TransformComponent object.
   *
   * @return The result.
   */
  Transform* duplicate() const;
  /**
   * @brief Copies the state of another TransformComponent into this object.
   *
   * @param another The TransformComponent instance to copy from.
   */
  void make_copy(const Transform &another);
  /**
   * @brief Assignment operator for TransformComponent.
   *
   * @param another The TransformComponent instance to copy from.
   */
  void operator=(const Transform &another);
  
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
  
protected: // Things that can be accessed in this class and subclasses.
  
  // stored here
  /** @brief The components. */
  std::vector<TransformComponent*> components;
  
  // points to one of the components
  /** @brief The root. */
  TransformComponent* root;
  
  // points to one of the components
  /** @brief The leaf. */
  TransformComponent* leaf;
  
};
}

#endif
