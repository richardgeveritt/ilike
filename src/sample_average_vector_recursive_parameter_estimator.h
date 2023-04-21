#ifndef SAMPLEAVERAGEVECTORRECURSIVEPARAMETERESTIMATOR_H
#define SAMPLEAVERAGEVECTORRECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

//#include <vector>
#include "vector_recursive_parameter_estimator.h"
#include "ilike_header.h"
#include "particle.h"

class SampleAverageVectorRecursiveParameterEstimator : public VectorRecursiveParameterEstimator
{

public:

  SampleAverageVectorRecursiveParameterEstimator();

  virtual ~SampleAverageVectorRecursiveParameterEstimator();

  SampleAverageVectorRecursiveParameterEstimator(const SampleAverageVectorRecursiveParameterEstimator &another);

  void operator=(const SampleAverageVectorRecursiveParameterEstimator &another);
  RecursiveParameterEstimator* duplicate() const;
  VectorRecursiveParameterEstimator* vector_duplicate() const;
  
  void update(const std::string &variable_name,
              const Particle &latest_particle,
              size_t iteration_counter,
              ProposalKernel* proposal);
  
protected:
  
  GainPtr gain;

  void make_copy(const SampleAverageVectorRecursiveParameterEstimator &another);

};

#endif
