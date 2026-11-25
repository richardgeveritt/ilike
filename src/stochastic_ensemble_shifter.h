#ifndef STOCHASTICENSEMBLESHIFTER_H
#define STOCHASTICENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{

class StochasticEnsembleShifter : public EnsembleShifter
{
  
public:
  
  StochasticEnsembleShifter();
  
  virtual ~StochasticEnsembleShifter();
  
  StochasticEnsembleShifter(const StochasticEnsembleShifter &another);
  
  void operator=(const StochasticEnsembleShifter &another);
  EnsembleShifter* duplicate() const;
  
  void setup(Ensemble* ensemble,
             double inverse_incremental_temperature);
  
  bool uses_shift_terms() const { return true; }
  void shift(const std::vector<arma::colvec> &shift_terms,
             arma::colvec &position,
             const std::vector<arma::colvec*> &measurements,
             const std::vector<arma::mat> &kalman_gains,
             const std::vector<arma::colvec> &myys,
             double inverse_incremental_temperature) const;
  
protected:
  
  void make_copy(const StochasticEnsembleShifter &another);
  
};
}

#endif
