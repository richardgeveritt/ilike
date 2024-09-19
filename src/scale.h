#ifndef SCALE_H
#define SCALE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

namespace ilike
{
class Scale
{
public:
  
  Scale();
  ~Scale();
  
  Scale(double constant_in);
  
  Scale(double constant_in,
        size_t dimension_in);
  
  Scale(const Scale &another);
  Scale& operator=(const Scale &another);
  
  Scale(Scale &&another);
  Scale& operator=(Scale &&another);
  
  double operator()() const;
  
  double get_constant() const;
  double& get_constant();
  
  //double operator()(size_t dimension);
  
private:
  
  void make_copy(const Scale &another);
  void make_copy(Scale &&another);
  
  double constant;
  bool divide_by_dimension;
  double dimension;
};
}

#endif
