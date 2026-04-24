#ifndef IIDDATASUBSETTER_H
#define IIDDATASUBSETTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "data_subsetter.h"

namespace ilike
{
  /**
   * @file iid_data_subsetter.h
   * @brief Defines the IIDDataSubsetter class.
   *
   * Implements iid data subsetter functionality derived from DataSubsetter. See the base class documentation for the full interface.
   *
   * @namespace ilike
   * @class IIDDataSubsetter
   * @brief An iid data subsetter derived from DataSubsetter.
   */


class IIDDataSubsetter : public DataSubsetter
{
  
public:
  
  /**
   * @brief Default constructor for IIDDataSubsetter.
   */
  IIDDataSubsetter();
  
  /**
   * @brief Destructor for IIDDataSubsetter.
   */
  virtual ~IIDDataSubsetter();
  
  /**
   * @brief Copy constructor for IIDDataSubsetter.
   *
   * @param another The IIDDataSubsetter instance to copy from.
   */
  IIDDataSubsetter(const IIDDataSubsetter &another);
  
  /**
   * @brief Assignment operator for IIDDataSubsetter.
   *
   * @param another The IIDDataSubsetter instance to copy from.
   */
  void operator=(const IIDDataSubsetter &another);
  /**
   * @brief Creates a deep copy of this IIDDataSubsetter object.
   *
   * @return The result.
   */
  DataSubsetter* duplicate() const;
  
  /**
   * @brief Performs the subset operation.
   *
   * @param num_pieces The num pieces.
   */
  void subset(size_t num_pieces);
  
protected:
  
  /**
   * @brief Copies the state of another IIDDataSubsetter into this object.
   *
   * @param another The IIDDataSubsetter instance to copy from.
   */
  void make_copy(const IIDDataSubsetter &another);
  
};
}

#endif
