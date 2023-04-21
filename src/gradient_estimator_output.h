#ifndef GRADIENTESTIMATOROUTPUT_H
#define GRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include "particle.h"
#include "parameters.h"

class GradientEstimatorOutput
{

public:

  GradientEstimatorOutput();
  virtual ~GradientEstimatorOutput();

  GradientEstimatorOutput(const GradientEstimatorOutput &another);

  void operator=(const GradientEstimatorOutput &another);
  virtual GradientEstimatorOutput* duplicate() const=0;
  
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Index* index,
                                        Particle &particle)=0;
  
  /*
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Index* index,
                                        Particle &particle,
                                        const Parameters &conditioned_on_parameters_in)=0;
  */
  
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Index* index,
                                                  Particle &particle)=0;
  
  /*
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Index* index,
                                                  Particle &particle,
                                                  const Parameters &conditioned_on_parameters_in)=0;
  */

protected:

  void make_copy(const GradientEstimatorOutput &another);

};

#endif
