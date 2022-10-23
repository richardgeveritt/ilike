#ifndef SCALE_H
#define SCALE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "particle.h"
#include "distributions.h"

class Scale
{
public:
  
  Scale();
  ~Scale();
  
  Scale(double constant_in);
  
  Scale(double constant_in,
        size_t dimension_in);
  
  Scale(const Scale &another);
  
  void operator=(const Scale &another);
  
  double operator()() const;
  
  double get_constant() const;
  double& get_constant();
  
  //double operator()(size_t dimension);
  
private:
  
  void make_copy(const Scale &another);
  
  double constant;
  bool divide_by_dimension;
  double dimension;
};

#endif
