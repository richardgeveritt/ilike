#ifndef GAUSSIANABCKERNELFACTOR_H
#define GAUSSIANABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "abc_kernel_factor.h"

namespace ilike
{
namespace exact_factor
{

class GaussianABCKernelFactor : public ABCKernelFactor
{
  
public:
  
  GaussianABCKernelFactor();
  GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                          const std::string &epsilon_variable_in,
                          Data* data_in);
  GaussianABCKernelFactor(const std::vector<std::string> &data_variables_in,
                          const std::string &epsilon_variable_in,
                          const std::string &weight_variable_in,
                          Data* data_in);
  
  virtual ~GaussianABCKernelFactor();
  
  GaussianABCKernelFactor(const GaussianABCKernelFactor &another);
  
  void operator=(const GaussianABCKernelFactor &another);
  Factor* duplicate() const;
  LikelihoodFactor* likelihood_factor_duplicate() const;
  ABCKernelFactor* abc_kernel_factor_duplicate() const;
  
  double likelihood_evaluate(const Parameters &input) const;
  
  void find_distance(const Parameters &input,
                     double &distance,
                     double &scale_constant) const;
  
  double evaluate_kernel_given_distance(const Parameters &input,
                                        double distance,
                                        double scale_constant) const;
  
  arma::mat likelihood_evaluate_gradient(const std::string &variable,
                                         const Parameters &input) const;
  
protected:
  
  void specific_set_data();
  
  void make_copy(const GaussianABCKernelFactor &another);
  
  double constant;
  
};
}
}

#endif
