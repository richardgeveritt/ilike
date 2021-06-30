#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"
#include "particles.h"

class SMCOutput : public LikelihoodEstimatorOutput
{

public:

  SMCOutput();
  virtual ~SMCOutput();

  SMCOutput(const SMCOutput &another);
  void operator=(const SMCOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

protected:

  std::deque<Particles> all_particles;

  std::deque<Particles> all_proposed;

  std::deque<arma::colvec> log_normalised_weights;

  std::deque<double> log_normalising_constant;

  void make_copy(const SMCOutput &another);

};

#endif
