#ifndef HMMINDEX_H
#define HMMINDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "index.h"

class HMMIndex : public Index
{

public:

  HMMIndex();
  HMMIndex(const std::vector<size_t> &factor_indices_time_zero_in,
           //bool evaluate_prior_in,
           const std::vector<size_t> &factor_indices_time_not_zero_in,
           bool evaluate_transition_model_in,
           size_t first_time_index_in);
  //HMMIndex(size_t single_index_in);

  virtual ~HMMIndex();

  HMMIndex(const HMMIndex &another);

  void operator=(const HMMIndex &another);
  Index* duplicate() const;
  HMMIndex* hmm_index_duplicate() const;
  
  std::vector<size_t>::const_iterator begin() const;
  std::vector<size_t>::const_iterator end() const;
  
  size_t size() const;
  
  arma::uvec get_uvec() const;
  
  ///std::vector<size_t> get_prior_indices() const;
  //std::vector<size_t> get_transition_indices() const;
  //std::vector<size_t> get_likelihood_indices() const;
  
  void set_time_index(size_t time_index_in);
  
  bool get_transition_model() const;
  
  //void add_index(const size_t &number);

protected:

  void make_copy(const HMMIndex &another);
  
  size_t time_index;
  size_t first_time_index;
  
  // If time = 0, which parts of HMMFactors are we evaluating?
  std::vector<size_t> likelihood_indices_time_zero;
  //bool evaluate_prior;
  
  // If time is not zero, which parts of HMMFactors are we evaluating?
  std::vector<size_t> likelihood_indices_time_not_zero;
  bool evaluate_transition_model;

};

#endif
