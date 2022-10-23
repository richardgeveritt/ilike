#ifndef SINGLEINDEX_H
#define SINGLEINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "index.h"

class SingleIndex : public Index
{

public:

  SingleIndex();
  virtual ~SingleIndex();

  SingleIndex(const SingleIndex &another);

  void operator=(const SingleIndex &another);
  virtual SingleIndex* single_index_duplicate() const=0;

protected:

  void make_copy(const SingleIndex &another);
  
};

#endif
