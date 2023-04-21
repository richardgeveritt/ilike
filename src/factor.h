#ifndef FACTOR_H
#define FACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

class Factor
{

public:

  Factor();
  virtual ~Factor();

  Factor(const Factor &another);

  void operator=(const Factor &another);
  virtual Factor* duplicate() const=0;
  
  virtual double evaluate(const Parameters &input) const=0;
  virtual arma::mat evaluate_gradient(const std::string &variable,
                                      const Parameters &input) const=0;

protected:
  
  void make_copy(const Factor &another);

};

#endif
