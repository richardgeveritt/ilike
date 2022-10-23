#ifndef VECTORSINGLEINDEX_H
#define VECTORSINGLEINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "single_index.h"

class VectorSingleIndex : public SingleIndex
{

public:

  VectorSingleIndex();
  VectorSingleIndex(const std::vector<size_t> &indices_in);
  VectorSingleIndex(size_t start_in,
                    size_t end_in);
  VectorSingleIndex(size_t single_index_in);

  virtual ~VectorSingleIndex();

  VectorSingleIndex(const VectorSingleIndex &another);

  void operator=(const VectorSingleIndex &another);
  Index* duplicate() const;
  SingleIndex* single_index_duplicate() const;
  
  std::vector<size_t>::const_iterator begin() const;
  std::vector<size_t>::const_iterator end() const;
  
  size_t size() const;
  
  arma::uvec get_uvec() const;

protected:

  void make_copy(const VectorSingleIndex &another);
  
  std::vector<size_t> indices;

};

#endif
