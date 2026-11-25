#ifndef SMCCRITERION_H
#define SMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "particles.h"
#include "ensemble.h"

namespace ilike
{

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
  
  virtual void subsample_find_desired_criterion(SMCOutput* current_state)=0;

  virtual double operator()(const Ensemble &particles) const=0;
  
  virtual void find_desired_criterion(EnsembleKalmanOutput* current_state)=0;

  virtual bool always_positive() const=0;

  // (dev/plan.md 3.21 (c2''), F2) a threshold in particles, for an SMC whose
  // number of particles changed
  void scale_desired_criterion(double factor) { this->desired_criterion *= factor; }
  
protected:
  
  double desired_criterion;
  
  void make_copy(const SMCCriterion &another);
  
};
}

#endif
