#ifndef INDEX_H
#define INDEX_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace ilike
{
class Index
{
  
public:
  
  Index();
  virtual ~Index();
  
  Index(const Index &another);
  
  void operator=(const Index &another);
  virtual Index* duplicate() const=0;
  
  // ok for now: might want to change later
  virtual std::vector<size_t>::const_iterator begin() const=0;
  virtual std::vector<size_t>::const_iterator end() const=0;
  
  virtual size_t size() const=0;
  
  virtual arma::uvec get_uvec() const=0;
  
  virtual bool get_transition_model() const=0;
  
protected:
  
  void make_copy(const Index &another);
  
};
}

#endif
