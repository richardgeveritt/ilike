#ifndef ESSSMCCRITERION_H
#define ESSSMCCRITERION_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "smc_criterion.h"

class SMCOutput;
class EnsembleKalmanOutput;

class ESSSMCCriterion : public SMCCriterion
{

public:

  ESSSMCCriterion();
  
  ESSSMCCriterion(double desired_criterion_in);

  virtual ~ESSSMCCriterion();

  ESSSMCCriterion(const ESSSMCCriterion &another);

  void operator=(const ESSSMCCriterion &another);
  SMCCriterion* duplicate() const;

  double operator()(const Particles &particles) const;
  double operator()(const Ensemble &particles) const;
  
  void find_desired_criterion(SMCOutput* current_state);
  
  void subsample_find_desired_criterion(SMCOutput* current_state);
  
  /*
  void find_desired_criterion(SMCOutput* current_state,
                              const Parameters &conditioned_on_parameters);
  void subsample_find_desired_criterion(SMCOutput* current_state,
                                        const Parameters &conditioned_on_parameters);
  */
  
  void find_desired_criterion(EnsembleKalmanOutput* current_state);
  
  /*
  void find_desired_criterion(EnsembleKalmanOutput* current_state,
                              const Parameters &conditioned_on_parameters);
  */
  
  bool always_positive() const;
  
protected:

  void make_copy(const ESSSMCCriterion &another);

};

#endif
