#ifndef ENSEMBLESHIFTER_H
#define ENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
//#include "ensemble_member.h"

namespace ilike
{
class Ensemble;
class EnsembleFactorVariables;

class EnsembleShifter
{
  
public:
  
  EnsembleShifter();
  virtual ~EnsembleShifter();
  
  EnsembleShifter(const EnsembleShifter &another);
  
  void operator=(const EnsembleShifter &another);
  virtual EnsembleShifter* duplicate() const=0;
  
  virtual void setup(Ensemble* ensemble,
                     double inverse_incremental_temperature)=0;
  
  virtual void shift(const EnsembleFactorVariables* ensemble_factor_variables,
                     arma::colvec &position,
                     const std::vector<arma::colvec*> &measurements,
                     const std::vector<arma::mat> &kalman_gains,
                     const std::vector<arma::colvec> &myys,
                     double inverse_incremental_temperature) const=0;
  
protected:
  
  void make_copy(const EnsembleShifter &another);
  
};
}

#endif
