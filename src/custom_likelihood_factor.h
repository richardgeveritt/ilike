#ifndef CUSTOMLIKELIHOODFACTOR_H
#define CUSTOMLIKELIHOODFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "ilike_header.h"
#include "likelihood_factor.h"

namespace ilike
{

class CustomLikelihoodFactor : public LikelihoodFactor
{
  
public:
  
  CustomLikelihoodFactor();
  CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                         Data* data_in);
  CustomLikelihoodFactor(EvaluateLogLikelihoodPtr likelihood_in,
                         EvaluateGradientLogLikelihoodPtr likelihood_gradient_in,
                         Data* data_in);
  
  virtual ~CustomLikelihoodFactor();
  
  CustomLikelihoodFactor(const CustomLikelihoodFactor &another);
  
  void operator=(const CustomLikelihoodFactor &another);
  Factor* duplicate() const;
  LikelihoodFactor* likelihood_factor_duplicate() const;
  
  double likelihood_evaluate(const Parameters &input) const;
  
  void set_maximum(EvaluateLogLikelihoodPtr maximum_in) { this->maximum = maximum_in; }
  double log_maximum(const Parameters &input) const
  {
    return (this->maximum != NULL) ? this->maximum(input, *this->data) : arma::datum::inf;
  }
  arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                         const Parameters &input) const;
  
protected:
  
  void specific_set_data();
  
  void make_copy(const CustomLikelihoodFactor &another);
  
  EvaluateLogLikelihoodPtr likelihood;
  EvaluateGradientLogLikelihoodPtr likelihood_gradient;
  EvaluateLogLikelihoodPtr maximum = NULL; // (density(…, maximum = cpp::f), dev/plan.md 3.21)
  
};
}

#endif
