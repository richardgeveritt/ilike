#ifndef ExactLikelihoodEstimatorOutput_H
#define ExactLikelihoodEstimatorOutput_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"
#include "particles.h"

class ExactLikelihoodEstimator;

class ExactLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{

public:

  ExactLikelihoodEstimatorOutput();
  ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in);
  virtual ~ExactLikelihoodEstimatorOutput();

  ExactLikelihoodEstimatorOutput(const ExactLikelihoodEstimatorOutput &another);
  void operator=(const ExactLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void continue_simulate(const Parameters &parameters);
  void estimate(const Parameters &parameters);

protected:

  // Stored in ModelAndAlgorithm.
  ExactLikelihoodEstimator* estimator;

  void make_copy(const ExactLikelihoodEstimatorOutput &another);

};

#endif
