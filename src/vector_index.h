#ifndef VECTORINDEX_H
#define VECTORINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "index.h"

namespace ilike
{
  /**
   * @file vector_index.h
   * @brief Defines the VectorIndex class.
   *
   * Implements vector index functionality derived from Index. See the base class documentation for the full interface.
   *
   * @namespace ilike
   * @class VectorIndex
   * @brief A vector index derived from Index.
   */


class VectorIndex : public Index
{
  
public:
  
  /**
   * @brief Default constructor for VectorIndex.
   */
  VectorIndex();
  /**
   * @brief Constructs a VectorIndex object.
   *
   * @param indices_in The indices.
   */
  VectorIndex(const std::vector<size_t> &indices_in);
  VectorIndex(size_t start_in,
              size_t end_in);
  /**
   * @brief Constructs a VectorIndex object.
   *
   * @param single_index_in The single index.
   */
  VectorIndex(size_t single_index_in);
  
  /**
   * @brief Destructor for VectorIndex.
   */
  virtual ~VectorIndex();
  
  /**
   * @brief Copy constructor for VectorIndex.
   *
   * @param another The VectorIndex instance to copy from.
   */
  VectorIndex(const VectorIndex &another);
  
  /**
   * @brief Assignment operator for VectorIndex.
   *
   * @param another The VectorIndex instance to copy from.
   */
  void operator=(const VectorIndex &another);
  /**
   * @brief Creates a deep copy of this VectorIndex object.
   *
   * @return The result.
   */
  Index* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a single_index pointer.
   *
   * @return The result.
   */
  Index* single_index_duplicate() const;
  
  /**
   * @brief Returns the indices.
   *
   * @return The result.
   */
  std::vector<size_t> get_indices() const;
  
  /**
   * @brief Returns the transition model.
   *
   * @return The result.
   */
  bool get_transition_model() const;
  
  /**
   * @brief Performs the begin operation.
   *
   * @return The result.
   */
  std::vector<size_t>::const_iterator begin() const;
  /**
   * @brief Performs the end operation.
   *
   * @return The result.
   */
  std::vector<size_t>::const_iterator end() const;
  
  /**
   * @brief Performs the size operation.
   *
   * @return The result.
   */
  size_t size() const;
  
  /**
   * @brief Returns the uvec.
   *
   * @return The result.
   */
  arma::uvec get_uvec() const;
  
  //void add_index(const size_t &number);
  
protected:
  
  /**
   * @brief Copies the state of another VectorIndex into this object.
   *
   * @param another The VectorIndex instance to copy from.
   */
  void make_copy(const VectorIndex &another);
  
  /** @brief The indices. */
  std::vector<size_t> indices;
  
};
}

#endif
