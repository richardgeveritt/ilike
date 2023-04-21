#ifndef ABCKERNELFACTOR_H
#define ABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "likelihood_factor.h"
#include "packing_instructions_map.h"

class ABCKernelFactor : public LikelihoodFactor
{

public:

  ABCKernelFactor();
  ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                  const std::string &epsilon_variable_in,
                  Data* data_in);
  ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                  const std::string &epsilon_variable_in,
                  const std::string &scale_variable_in,
                  Data* data_in);

  virtual ~ABCKernelFactor();

  ABCKernelFactor(const ABCKernelFactor &another);

  void operator=(const ABCKernelFactor &another);
  virtual ABCKernelFactor* abc_kernel_factor_duplicate() const=0;
  
  virtual void find_distance(const Parameters &input,
                             double &distance,
                             double &scale_constant) const=0;
  
  virtual double evaluate_kernel_given_distance(const Parameters &input,
                                                double distance,
                                                double scale_constant) const=0;

protected:
  
  void make_copy(const ABCKernelFactor &another);
  
  std::vector<std::string> data_variables;
  
  std::string epsilon_variable;
  
  std::string scale_variable;
  
  arma::colvec data_colvec;
  
  PackingInstructionsMap packing_instructions;

};

#endif
