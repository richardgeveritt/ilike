#ifndef STANDARDMCMCOUTPUT_H
#define STANDARDMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"
#include "parameters.h"

class StandardMCMCOutput : public MoveOutput
{

public:

  StandardMCMCOutput(const Parameters &parameter_in);
  StandardMCMCOutput();

  virtual ~StandardMCMCOutput();

  StandardMCMCOutput(const StandardMCMCOutput &another);

  void operator=(const StandardMCMCOutput &another);
  MoveOutput* duplicate() const;
  
  void push_back(const Particle &particle_in);

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
protected:

  void make_copy(const StandardMCMCOutput &another);
  
  std::deque<Particle> output;

};

#endif
