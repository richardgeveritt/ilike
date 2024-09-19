#ifndef ABCKERNELFACTOR_H
#define ABCKERNELFACTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "likelihood_factor.h"
#include "packing_instructions_map.h"

namespace ilike
{

/// An ABCKernel type of Factor, to be used in an ExactLikelihoodEstimator.
/// ExactLikelihoodEstimator contains LikelihoodFactor%s, DistributionFactor%s and IndependentProposalKernel%s. These will be evaluated
/// when a Particle is evaluating its factors in a weight update or Metropolis step, for example. The ABCKernelFactor has different subtypes,
/// corresponding to different types of ABC kernel. Each ABC kernel is a function of the distance between the observed and simulated data.
class ABCKernelFactor : public LikelihoodFactor
{
  
public:
  
  /// Default constructor.
  ABCKernelFactor();
  
  /// Constructor reading in the variables that we look at in the data, the tolerance epsilon, and the data.
  ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                  const std::string &epsilon_variable_in,
                  Data* data_in);
  
  /// Constructor reading in the variables that we look at in the data, the tolerance epsilon, the variable corresponding to the kernel scale and the data.
  ABCKernelFactor(const std::vector<std::string> &data_variables_in,
                  const std::string &epsilon_variable_in,
                  const std::string &scale_variable_in,
                  Data* data_in);
  
  /// Default destructor.
  virtual ~ABCKernelFactor();
  
  /// Copy constructor..
  ABCKernelFactor(const ABCKernelFactor &another);
  
  /// Equals operator.
  ABCKernelFactor& operator=(const ABCKernelFactor &another);
  
  /// Deep copy.
  virtual ABCKernelFactor* abc_kernel_factor_duplicate() const=0;
  
  /// Finds the distance and scaling for the kernel from the given input, returned by reference.
  virtual void find_distance(const Parameters &input,
                             double &distance,
                             double &scale_constant) const=0;
  
  /// Evaluates the kernel given the distance and scaling.
  virtual double evaluate_kernel_given_distance(const Parameters &input,
                                                double distance,
                                                double scale_constant) const=0;
  
protected:
  
  /// Does the work of copying the class.
  void make_copy(const ABCKernelFactor &another);
  
  /// The variables in the data used in the ABC kernel.
  std::vector<std::string> data_variables;
  
  /// The tolerance variable.
  std::string epsilon_variable;
  
  /// The scale variable.
  std::string scale_variable;
  
  /// The observed data as a column vector.
  arma::colvec data_colvec;
  
  /// The instructions as to how to pack the data variables into a column vector.
  PackingInstructionsMap packing_instructions;
  
};
}

#endif
