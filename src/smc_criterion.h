#ifndef SMCCRITERION_H
#define SMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"
#include "ensemble.h"

class SMCOutput;
class EnsembleKalmanOutput;

class SMCCriterion
{

public:

  SMCCriterion();
  virtual ~SMCCriterion();
  
  SMCCriterion(double desired_criterion_in);

  SMCCriterion(const SMCCriterion &another);

  void operator=(const SMCCriterion &another);
  virtual SMCCriterion* duplicate() const=0;

  virtual double operator()(const Particles &particles) const=0;
  
  virtual void find_desired_criterion(SMCOutput* current_state)=0;
  virtual void find_desired_criterion(SMCOutput* current_state,
                                      const Parameters &conditioned_on_parameters)=0;
  
  virtual double operator()(const Ensemble &particles) const=0;
  
  virtual void find_desired_criterion(EnsembleKalmanOutput* current_state)=0;
  virtual void find_desired_criterion(EnsembleKalmanOutput* current_state,
                                      const Parameters &conditioned_on_parameters)=0;

protected:

  double desired_criterion;

  void make_copy(const SMCCriterion &another);

};

#endif
