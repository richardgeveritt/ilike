#ifndef ENSEMBLESHIFTER_H
#define ENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

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
  
  // (C279) shift_terms: each active measurement factor's predicted
  // measurement plus its noise scaled by sqrt(inverse incremental
  // temperature), for the shifters that use them (uses_shift_terms); empty
  // otherwise. The caller computes them from the member's measurement outputs
  // or from its rows of a compact ensemble.
  virtual bool uses_shift_terms() const { return false; }
  virtual void shift(const std::vector<arma::colvec> &shift_terms,
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
