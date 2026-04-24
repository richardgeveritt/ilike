#ifndef DATASUBSAMPLER_H
#define DATASUBSAMPLER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>
#include "parameters.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file data_subsampler.h
   * @brief Defines the DataSubsetter class.
   *
   * Provides data subsetter functionality.
   *
   * @namespace ilike
   * @class DataSubsetter
   * @brief The data subsetter class.
   */



class DataSubsetter;
class IIDDataSubsetter;

class DataSubsampler
{
  
public:
  
  /**
   * @brief Performs the datasubsampler operation.
   */
  DataSubsampler();
  /**
   * @brief Performs the ~datasubsampler operation.
   */
  virtual ~DataSubsampler();
  
  /**
   * @brief Performs the datasubsampler operation.
   *
   * @param another The DataSubsetter instance to copy from.
   */
  DataSubsampler(const DataSubsampler &another);
  
  /**
   * @brief Assignment operator for DataSubsetter.
   *
   * @param another The DataSubsetter instance to copy from.
   */
  void operator=(const DataSubsampler &another);
  /**
   * @brief Creates a deep copy of this DataSubsetter object.
   *
   * @return The result.
   */
  virtual DataSubsampler* duplicate() const=0;
  
  /**
   * @brief Subsamples the data for use in the algorithm.
   *
   * @param num_pieces The num pieces.
   */
  void subsample(size_t num_pieces);
  
  // not implemented at the moment
  // need a boost::multi_index structure to replace the vector below
  //void subset(size_t num_pieces
  //            const std::string &variable);
  
  // Stored here.
  // At the moment, this is simply the case where we make a new dataset, copying elements of the actual data.
  // At some point, we want to subclass data, where we simply point to elements of the data.
  // Need to check if we can do this, and pass this class structure to R.
  // Need to check if we can still have the user write functions that take (const Data &data_in)
  Data* small_data;
  
  double ratio;
  
protected:
  
  friend IIDDataSubsetter;
  // Not stored here. Stored in "main'.
  /** @brief The rng. */
  RandomNumberGenerator* rng;
  
  /**
   * @brief Copies the state of another DataSubsetter into this object.
   *
   * @param another The DataSubsetter instance to copy from.
   */
  void make_copy(const DataSubsampler &another);
  
  /** @brief The variables. */
  std::vector<std::string> variables;
  
  // Stored here.
  /** @brief The subsetters. */
  std::vector<DataSubsetter*> subsetters;
  
};
}

#endif
