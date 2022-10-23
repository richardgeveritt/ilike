#ifndef PDMPMCMCOUTPUT_H
#define PDMPMCMCOUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "move_output.h"

class PDMPMCMCOutput : public MoveOutput
{

public:

  PDMPMCMCOutput();

  virtual ~PDMPMCMCOutput();

  PDMPMCMCOutput(const PDMPMCMCOutput &another);

  void operator=(const PDMPMCMCOutput &another);
  MoveOutput* duplicate() const;

  Particle& back();
  Particle back() const;
  
  std::vector<Parameters> get_vector_of_parameters() const;
  
protected:
  
  Particle dummy;

  void make_copy(const PDMPMCMCOutput &another);

};

#endif
