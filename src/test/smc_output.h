//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "particles.h"

#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

class SMCOutput
{

public:

  SMCOutput();

  virtual ~SMCOutput();

protected:

  std::deque<Particles> all_particles;

  std::deque<Particles> all_proposed;

  std::deque<arma::colvec> log_normalised_weights;

  std::deque<double> log_normalising_constant;

};

#endif
