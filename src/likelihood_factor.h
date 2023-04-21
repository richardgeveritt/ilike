#ifndef LIKELIHOODFACTOR_H
#define LIKELIHOODFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "factor.h"

class LikelihoodFactor : public Factor
{

public:

  LikelihoodFactor();
  LikelihoodFactor(Data* data_in);
  virtual ~LikelihoodFactor();

  LikelihoodFactor(const LikelihoodFactor &another);

  void operator=(const LikelihoodFactor &another);
  virtual LikelihoodFactor* likelihood_factor_duplicate() const=0;
  
  double evaluate(const Parameters &input) const;
  virtual double likelihood_evaluate(const Parameters &input) const=0;
  
  arma::mat evaluate_gradient(const std::string &variable,
                              const Parameters &input) const;
  virtual arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                                 const Parameters &input) const=0;
  
  void set_data(Data* data_in);

protected:
  
  virtual void specific_set_data()=0;

  void make_copy(const LikelihoodFactor &another);
  
  // Not stored here. Stored in "main'.
  Data* data;
  
};

#endif
