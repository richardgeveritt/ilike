#ifndef SINGLEPOINTMOVEOUTPUT_H
#define SINGLEPOINTMOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"
#include "particle.h"

class SinglePointMoveOutput : public MoveOutput
{

public:

  SinglePointMoveOutput(const Parameters &parameter_in);
  SinglePointMoveOutput(const Particle &particle_in);
  SinglePointMoveOutput();

  virtual ~SinglePointMoveOutput();

  SinglePointMoveOutput(const SinglePointMoveOutput &another);

  void operator=(const SinglePointMoveOutput &another);
  MoveOutput* duplicate() const;

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
protected:

  void make_copy(const SinglePointMoveOutput &another);
  
  Particle output;

};

#endif
