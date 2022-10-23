#ifndef MOVEOUTPUT_H
#define MOVEOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particle.h"

class MoveOutput
{

public:

  MoveOutput();
  virtual ~MoveOutput();

  MoveOutput(const MoveOutput &another);

  void operator=(const MoveOutput &another);
  virtual MoveOutput* duplicate() const=0;

  virtual Particle& back()=0;
  virtual Particle back() const=0;
  
  virtual std::vector<Parameters> get_vector_of_parameters() const=0;

protected:

  void make_copy(const MoveOutput &another);

};

#endif
