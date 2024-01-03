#ifndef VECTORINDEX_H
#define VECTORINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "index.h"

class VectorIndex : public Index
{

public:

  VectorIndex();
  VectorIndex(const std::vector<size_t> &indices_in);
  VectorIndex(size_t start_in,
                    size_t end_in);
  VectorIndex(size_t single_index_in);

  virtual ~VectorIndex();

  VectorIndex(const VectorIndex &another);

  void operator=(const VectorIndex &another);
  Index* duplicate() const;
  Index* single_index_duplicate() const;
  
  std::vector<size_t> get_indices() const;
  
  bool get_transition_model() const;
  
  std::vector<size_t>::const_iterator begin() const;
  std::vector<size_t>::const_iterator end() const;
  
  size_t size() const;
  
  arma::uvec get_uvec() const;
  
  //void add_index(const size_t &number);

protected:

  void make_copy(const VectorIndex &another);
  
  std::vector<size_t> indices;

};

#endif
