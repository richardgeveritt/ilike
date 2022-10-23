#ifndef SAMPLEAVERAGEGAUSSIANRECURSIVEPARAMETERESTIMATOR_H
#define SAMPLEAVERAGEGAUSSIANRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "gaussian_recursive_parameter_estimator.h"
#include "particle.h"
#include "function_pointers.h"

class ProposalKernel;

class SampleAverageGaussianRecursiveParameterEstimator : public GaussianRecursiveParameterEstimator
{

public:

  SampleAverageGaussianRecursiveParameterEstimator();

  virtual ~SampleAverageGaussianRecursiveParameterEstimator();

  SampleAverageGaussianRecursiveParameterEstimator(const SampleAverageGaussianRecursiveParameterEstimator &another);

  void operator=(const SampleAverageGaussianRecursiveParameterEstimator &another);
  RecursiveParameterEstimator* duplicate() const;
  GaussianRecursiveParameterEstimator* gaussian_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  GainPtr gain;

  void make_copy(const SampleAverageGaussianRecursiveParameterEstimator &another);

};

#endif
