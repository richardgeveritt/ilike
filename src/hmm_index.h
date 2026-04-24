#ifndef HMMINDEX_H
#define HMMINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "index.h"

namespace ilike
{
  /**
   * @file hmm_index.h
   * @brief Defines the HMMIndex class.
   *
   * Implements hmm index functionality derived from Index. See the base class documentation for the full interface.
   *
   * @namespace ilike
   * @class HMMIndex
   * @brief A hmm index derived from Index.
   */


class HMMIndex : public Index
{
  
public:
  
  /**
   * @brief Default constructor for HMMIndex.
   */
  HMMIndex();
  HMMIndex(const std::vector<size_t> &factor_indices_time_zero_in,
           //bool evaluate_prior_in,
           const std::vector<size_t> &factor_indices_time_not_zero_in,
           bool evaluate_transition_model_in,
           size_t first_time_index_in);
  //HMMIndex(size_t single_index_in);
  
  /**
   * @brief Destructor for HMMIndex.
   */
  virtual ~HMMIndex();
  
  /**
   * @brief Copy constructor for HMMIndex.
   *
   * @param another The HMMIndex instance to copy from.
   */
  HMMIndex(const HMMIndex &another);
  
  /**
   * @brief Assignment operator for HMMIndex.
   *
   * @param another The HMMIndex instance to copy from.
   */
  void operator=(const HMMIndex &another);
  /**
   * @brief Creates a deep copy of this HMMIndex object.
   *
   * @return The result.
   */
  Index* duplicate() const;
  /**
   * @brief Creates a deep copy and returns it as a hmm_index pointer.
   *
   * @return The result.
   */
  HMMIndex* hmm_index_duplicate() const;
  
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
  
  ///std::vector<size_t> get_prior_indices() const;
  //std::vector<size_t> get_transition_indices() const;
  //std::vector<size_t> get_likelihood_indices() const;
  
  /**
   * @brief Sets the time index.
   *
   * @param time_index_in The time index.
   */
  void set_time_index(size_t time_index_in);
  
  /**
   * @brief Returns the transition model.
   *
   * @return The result.
   */
  bool get_transition_model() const;
  
  //void add_index(const size_t &number);
  
protected:
  
  /**
   * @brief Copies the state of another HMMIndex into this object.
   *
   * @param another The HMMIndex instance to copy from.
   */
  void make_copy(const HMMIndex &another);
  
  /** @brief The time index. */
  size_t time_index;
  /** @brief The first time index. */
  size_t first_time_index;
  
  // If time = 0, which parts of HMMFactors are we evaluating?
  /** @brief The likelihood indices time zero. */
  std::vector<size_t> likelihood_indices_time_zero;
  //bool evaluate_prior;
  
  // If time is not zero, which parts of HMMFactors are we evaluating?
  /** @brief The likelihood indices time not zero. */
  std::vector<size_t> likelihood_indices_time_not_zero;
  /** @brief The evaluate transition model. */
  bool evaluate_transition_model;
  
};
}

#endif
