#include <RcppArmadillo.h>
using namespace Rcpp;

#ifndef SEQUENCER_H
#define SEQUENCER_H

class Sequencer
{
public:

  Sequencer(void);
  Sequencer(const Sequencer &another);
  virtual ~Sequencer(void);

  void operator=(const Sequencer &another);
  virtual Sequencer* duplicate() const=0;

protected:

  void make_copy(const Sequencer &another);

  unsigned int current_phase;
  double current_target;

  std::vector< std::vector<double> > schedule;

  std::vector<double> desired_cess;

};

#endif
