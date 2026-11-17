#ifndef KALMANPREDICTOR_H
#define KALMANPREDICTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{

class KalmanFilterOutput;

class KalmanPredictor
{
  
public:
  
  KalmanPredictor();
  virtual ~KalmanPredictor();
  
  KalmanPredictor(const KalmanPredictor &another);
  
  void operator=(const KalmanPredictor &another);
  virtual KalmanPredictor* duplicate() const=0;
  
  virtual void predict(KalmanFilterOutput* current_state)=0;

  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;
  
protected:
  
  bool set_using_parameters;
  Parameters conditioned_on_parameters;
  
  void make_copy(const KalmanPredictor &another);
  
};
}

#endif
