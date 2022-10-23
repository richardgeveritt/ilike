#ifndef IIDDATASUBSETTER_H
#define IIDDATASUBSETTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "data_subsetter.h"

class IIDDataSubsetter : public DataSubsetter
{

public:

  IIDDataSubsetter();

  virtual ~IIDDataSubsetter();

  IIDDataSubsetter(const IIDDataSubsetter &another);

  void operator=(const IIDDataSubsetter &another);
  DataSubsetter* duplicate() const;
  
  void subset(size_t num_pieces);
  
protected:

  void make_copy(const IIDDataSubsetter &another);

};

#endif
