#ifndef INDEX_H
#define INDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace ilike
{
  /**
   * @file index.h
   * @brief Defines the Index class.
   *
   * Provides index functionality.
   *
   * @namespace ilike
   * @class Index
   * @brief The index class.
   */


class Index
{
  
public:
  
  /**
   * @brief Default constructor for Index.
   */
  Index();
  /**
   * @brief Destructor for Index.
   */
  virtual ~Index();
  
  /**
   * @brief Copy constructor for Index.
   *
   * @param another The Index instance to copy from.
   */
  Index(const Index &another);
  
  /**
   * @brief Assignment operator for Index.
   *
   * @param another The Index instance to copy from.
   */
  void operator=(const Index &another);
  /**
   * @brief Creates a deep copy of this Index object.
   *
   * @return The result.
   */
  virtual Index* duplicate() const=0;
  
  // ok for now: might want to change later
  /**
   * @brief Performs the begin operation.
   *
   * @return The result.
   */
  virtual std::vector<size_t>::const_iterator begin() const=0;
  /**
   * @brief Performs the end operation.
   *
   * @return The result.
   */
  virtual std::vector<size_t>::const_iterator end() const=0;
  
  /**
   * @brief Performs the size operation.
   *
   * @return The result.
   */
  virtual size_t size() const=0;
  
  /**
   * @brief Returns the uvec.
   *
   * @return The result.
   */
  virtual arma::uvec get_uvec() const=0;
  
  /**
   * @brief Returns the transition model.
   *
   * @return The result.
   */
  virtual bool get_transition_model() const=0;
  
protected:
  
  /**
   * @brief Copies the state of another Index into this object.
   *
   * @param another The Index instance to copy from.
   */
  void make_copy(const Index &another);
  
};
}

#endif
