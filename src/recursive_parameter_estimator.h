#ifndef RECURSIVEPARAMETERESTIMATOR_H
#define RECURSIVEPARAMETERESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particle.h"
class ProposalKernel;

class RecursiveParameterEstimator
{

public:

  RecursiveParameterEstimator();
  virtual ~RecursiveParameterEstimator();

  RecursiveParameterEstimator(const RecursiveParameterEstimator &another);

  void operator=(const RecursiveParameterEstimator &another);
  
  virtual void update(const std::string &variable_name,
                      const Particle &latest_particle,
                      size_t iteration_counter,
                      ProposalKernel* proposal)=0;

protected:

  void make_copy(const RecursiveParameterEstimator &another);

};

#endif
