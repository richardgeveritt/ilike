#include "square_root_ensemble_shifter.h"
#include "ensemble_factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble.h"

// follows https://arxiv.org/abs/2006.02941

namespace ilike
{
SquareRootEnsembleShifter::SquareRootEnsembleShifter()
:EnsembleShifter()
{
}

SquareRootEnsembleShifter::~SquareRootEnsembleShifter()
{
  
}

//Copy constructor for the SquareRootEnsembleShifter class.
SquareRootEnsembleShifter::SquareRootEnsembleShifter(const SquareRootEnsembleShifter &another)
:EnsembleShifter(another)
{
  this->make_copy(another);
}

void SquareRootEnsembleShifter::operator=(const SquareRootEnsembleShifter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  EnsembleShifter::operator=(another);
  this->make_copy(another);
}

EnsembleShifter* SquareRootEnsembleShifter::duplicate() const
{
  return( new SquareRootEnsembleShifter(*this));
}

void SquareRootEnsembleShifter::make_copy(const SquareRootEnsembleShifter &another)
{
  this->mean_position = another.mean_position;
  this->As = another.As;
}

void SquareRootEnsembleShifter::setup(Ensemble* ensemble,
                                      double inverse_incremental_temperature)
{
  this->mean_position = arma::mean(ensemble->packed_members,0).as_col();

  this->As = ensemble->ensemble_factors->get_sqrt_adjustments(ensemble->Cxys,
                                                              ensemble->Cyys,
                                                              inverse_incremental_temperature);
  
}

void SquareRootEnsembleShifter::shift(const EnsembleFactorVariables* ensemble_factor_variables,
                                      arma::colvec &position,
                                      const std::vector<arma::colvec*> &measurements,
                                      const std::vector<arma::mat> &kalman_gains,
                                      const std::vector<arma::colvec> &myys,
                                      double inverse_incremental_temperature) const
{
  std::vector<arma::colvec> shift_terms = ensemble_factor_variables->get_deterministic_shifts();

  arma::colvec second_part_shiftee = position - this->mean_position;
  
  // shift taken from https://arxiv.org/abs/2204.04386
  for (size_t j=0;
       j<kalman_gains.size();
       ++j)
  {
    
    position = this->mean_position + kalman_gains[j]*(*measurements[j] - myys[j]) + this->As[j]*second_part_shiftee;
  }
}
}
