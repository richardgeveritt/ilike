#ifndef DATASUBSETTER_H
#define DATASUBSETTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "distributions.h"
#include "data.h"

class DataSubsampler;

class DataSubsetter
{

public:

  DataSubsetter();
  virtual ~DataSubsetter();

  DataSubsetter(const DataSubsetter &another);

  void operator=(const DataSubsetter &another);
  virtual DataSubsetter* duplicate() const=0;
  
  virtual void subset(size_t num_pieces)=0;

protected:

  void make_copy(const DataSubsetter &another);
  
  // Not stored here. Stored in "main'.
  Data* data;
  
  // not stored here
  DataSubsampler* subsampler;

};

#endif
