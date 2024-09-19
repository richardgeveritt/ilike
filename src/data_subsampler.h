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

class DataSubsetter;
class IIDDataSubsetter;

class DataSubsampler
{
  
public:
  
  DataSubsampler();
  virtual ~DataSubsampler();
  
  DataSubsampler(const DataSubsampler &another);
  
  void operator=(const DataSubsampler &another);
  virtual DataSubsampler* duplicate() const=0;
  
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
  RandomNumberGenerator* rng;
  
  void make_copy(const DataSubsampler &another);
  
  std::vector<std::string> variables;
  
  // Stored here.
  std::vector<DataSubsetter*> subsetters;
  
};
}

#endif
