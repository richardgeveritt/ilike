#ifndef DATASUBSETTER_H
#define DATASUBSETTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distributions.h"
#include "parameters.h"

namespace ilike
{
  /**
   * @file data_subsetter.h
   * @brief Defines the DataSubsampler class.
   *
   * Provides data subsampler functionality.
   *
   * @namespace ilike
   * @class DataSubsampler
   * @brief The data subsampler class.
   */


class DataSubsampler;

class DataSubsetter
{
  
public:
  
  /**
   * @brief Performs the datasubsetter operation.
   */
  DataSubsetter();
  /**
   * @brief Performs the ~datasubsetter operation.
   */
  virtual ~DataSubsetter();
  
  /**
   * @brief Performs the datasubsetter operation.
   *
   * @param another The DataSubsampler instance to copy from.
   */
  DataSubsetter(const DataSubsetter &another);
  
  /**
   * @brief Assignment operator for DataSubsampler.
   *
   * @param another The DataSubsampler instance to copy from.
   */
  void operator=(const DataSubsetter &another);
  /**
   * @brief Creates a deep copy of this DataSubsampler object.
   *
   * @return The result.
   */
  virtual DataSubsetter* duplicate() const=0;
  
  /**
   * @brief Performs the subset operation.
   *
   * @param num_pieces The num pieces.
   */
  virtual void subset(size_t num_pieces)=0;
  
protected:
  
  /**
   * @brief Copies the state of another DataSubsampler into this object.
   *
   * @param another The DataSubsampler instance to copy from.
   */
  void make_copy(const DataSubsetter &another);
  
  // Not stored here. Stored in "main'.
  /** @brief The data. */
  Data* data;
  
  // not stored here
  /** @brief The subsampler. */
  DataSubsampler* subsampler;
  
};
}

#endif
