#ifndef PROPOSALSTORE_H
#define PROPOSALSTORE_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "parameters.h"

namespace ilike
{
class GradientEstimatorOutput;

class ProposalStore
{
  
public:
  
  ProposalStore();
  virtual ~ProposalStore();
  
  ProposalStore(const Parameters &transformed_parameters_in);
  
  ProposalStore(GradientEstimatorOutput* gradient_estimator_output_in);
  
  ProposalStore(const Parameters &transformed_parameters_in,
                GradientEstimatorOutput* gradient_estimator_output_in);
  
  ProposalStore(const ProposalStore &another);
  
  void operator=(const ProposalStore &another);
  
  void set_transformed_parameters(const Parameters &transformed_parameters_in);
  void set_gradient_estimator_output(GradientEstimatorOutput* gradient_estimator_output_in);
  
  GradientEstimatorOutput* get_gradient_estimator_output() const;
  
  Parameters get_transformed_parameters() const;
  
protected:
  
  Parameters transformed_parameters;
  
  // stored here
  GradientEstimatorOutput* gradient_estimator_output;
  
  void make_copy(const ProposalStore &another);
  
};
}

#endif
