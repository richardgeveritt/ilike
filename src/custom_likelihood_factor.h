#ifndef CUSTOMLIKELIHOODFACTOR_H
#define CUSTOMLIKELIHOODFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ilike_header.h"
#include "likelihood_factor.h"

class CustomLikelihoodFactor : public LikelihoodFactor
{

public:

  CustomLikelihoodFactor();
  CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                         Data* data_in);

  virtual ~CustomLikelihoodFactor();

  CustomLikelihoodFactor(const CustomLikelihoodFactor &another);

  void operator=(const CustomLikelihoodFactor &another);
  Factor* duplicate() const;
  LikelihoodFactor* likelihood_factor_duplicate() const;
  
  double likelihood_evaluate(const Parameters &input) const;
  
  arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                         const Parameters &input) const;

protected:
  
  void specific_set_data();

  void make_copy(const CustomLikelihoodFactor &another);
  
  EvaluateLogLikelihoodPtr likelihood;
  EvaluateGradientLogLikelihoodPtr likelihood_gradient;

};

#endif
