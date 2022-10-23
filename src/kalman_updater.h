#ifndef KALMANUPDATER_H
#define KALMANUPDATER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

class KalmanFilterOutput;
class KalmanFilter;

class KalmanUpdater
{

public:

  KalmanUpdater();
  virtual ~KalmanUpdater();

  KalmanUpdater(const KalmanUpdater &another);

  void operator=(const KalmanUpdater &another);
  virtual KalmanUpdater* duplicate() const=0;

  virtual void update(KalmanFilterOutput* current_state,
                      const arma::colvec &current_measurement)=0;
  
  virtual void set_parameters(const Parameters &conditioned_on_parameters_in)=0;

protected:
  
  friend KalmanFilter;
  bool set_using_parameters;
  
  void make_copy(const KalmanUpdater &another);

};

#endif
