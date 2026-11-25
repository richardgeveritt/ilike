#ifndef SQUAREROOTENSEMBLESHIFTER_H
#define SQUAREROOTENSEMBLESHIFTER_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "ensemble_shifter.h"

namespace ilike
{

class SquareRootEnsembleShifter : public EnsembleShifter
{
  
public:
  
  SquareRootEnsembleShifter();
  
  virtual ~SquareRootEnsembleShifter();
  
  SquareRootEnsembleShifter(const SquareRootEnsembleShifter &another);
  
  void operator=(const SquareRootEnsembleShifter &another);
  EnsembleShifter* duplicate() const;
  
  void setup(Ensemble* ensemble,
             double inverse_incremental_temperature);
  
  void shift(const std::vector<arma::colvec> &shift_terms,
             arma::colvec &position,
             const std::vector<arma::colvec*> &measurements,
             const std::vector<arma::mat> &kalman_gains,
             const std::vector<arma::colvec> &myys,
             double inverse_incremental_temperature) const;
  
protected:
  
  std::vector<arma::mat> As;
  arma::colvec mean_position;

  void make_copy(const SquareRootEnsembleShifter &another);
  
};
}

#endif
