#ifndef POSITIVESMCCRITERION_H
#define POSITIVESMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_criterion.h"

namespace ilike
{
class PositiveSMCCriterion : public SMCCriterion
{
  
public:
  
  PositiveSMCCriterion();
  
  virtual ~PositiveSMCCriterion();
  
  PositiveSMCCriterion(const PositiveSMCCriterion &another);
  
  void operator=(const PositiveSMCCriterion &another);
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
  
  void make_copy(const PositiveSMCCriterion &another);
  
};
}

#endif
