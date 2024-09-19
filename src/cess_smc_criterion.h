#ifndef CESSSMCCRITERION_H
#define CESSSMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_criterion.h"

namespace ilike
{
class SMCOutput;
class EnsembleKalmanOutput;

class CESSSMCCriterion : public SMCCriterion
{
  
public:
  
  CESSSMCCriterion();
  
  CESSSMCCriterion(double desired_criterion_in);
  
  virtual ~CESSSMCCriterion();
  
  CESSSMCCriterion(const CESSSMCCriterion &another);
  
  void operator=(const CESSSMCCriterion &another);
  SMCCriterion* duplicate() const;
  
  double operator()(const Particles &particles) const;
  
  void find_desired_criterion(SMCOutput* current_state);
  
  void subsample_find_desired_criterion(SMCOutput* current_state);
  /*
   void find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters);
   void subsample_find_desired_criterion(SMCOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  double operator()(const Ensemble &particles) const;
  
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  /*
   void find_desired_criterion(EnsembleKalmanOutput* current_state,
   const Parameters &conditioned_on_parameters);
   */
  
  bool always_positive() const;
  
protected:
  
  void make_copy(const CESSSMCCriterion &another);
  
};
}

#endif
