#include "square_root_ensemble_shifter.h"
#include "ensemble_factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble.h"

// follows https://arxiv.org/abs/2006.02941

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
  this->Sigma = another.Sigma;
  this->HSigmaHts = another.HSigmaHts;
  this->As = another.As;
}

void SquareRootEnsembleShifter::setup(Ensemble* ensemble,
                                      double inverse_incremental_temperature)
{
  this->mean_position = arma::mean(ensemble->packed_members,0).as_col();
  this->Sigma = arma::cov(ensemble->packed_members);
  
  this->HSigmaHts.clear();
  this->HSigmaHts.reserve(ensemble->Cxys.size());
  for (size_t i=0;
       i<ensemble->Cxys.size();
       ++i)
  {
    this->HSigmaHts.push_back(arma::cov(ensemble->packed_measurement_states[i]));
  }
  
  this->As = ensemble->ensemble_factors->get_sqrt_adjustments(this->Sigma,
                                                              this->HSigmaHts,
                                                              inverse_incremental_temperature);
  
  /*
  arma::mat U;
  arma::vec Gdiag;
  arma::svd(this->Ftranspose,Gdiag,U,this->Zf);
  this->Ftranspose = this->Ftranspose.t();
  arma::mat G(this->Zf.n_rows,this->Zf.n_cols);
  G.diag() = Gdiag;
  this->Ginv = arma::pinv(G);
  
  this->Vs.clear();
  this->Vs.reserve(ensemble->Cxys.size());
  for (size_t i=0;
       i<ensemble->Cxys.size();
       ++i)
  {
    this->Vs.push_back(ensemble->packed_measurement_states[i].each_row()-arma::conv_to<arma::rowvec>::from(ensemble->myys[i]));
  }
  
  As = ensemble->ensemble_factors->get_adjustments(this->Zf,
                                                   this->Ginv,
                                                   this->Ftranspose,
                                                   this->Vs,
                                                   inverse_incremental_temperature);
  */
}

void SquareRootEnsembleShifter::shift(const EnsembleFactorVariables* ensemble_factor_variables,
                                      arma::colvec &position,
                                      const std::vector<arma::colvec*> &measurements,
                                      const std::vector<arma::mat> &kalman_gains,
                                      const std::vector<arma::colvec> &myys,
                                      double inverse_incremental_temperature) const
{
  std::vector<arma::colvec> shift_terms = ensemble_factor_variables->get_deterministic_shifts();
  
  //std::vector<arma::colvec> first_part_of_shift;
  //first_part_of_shift.reserve(Cxys.size());
  
  arma::colvec second_part_shiftee = position - this->mean_position;
  //position = this->mean_position;
  
  // shift taken from https://arxiv.org/abs/2204.04386
  for (size_t j=0;
       j<kalman_gains.size();
       ++j)
  {
    //position = position + kalman_gains[j]*(*measurements[j] - shift_terms[j]) + this->As[j]*second_part_shiftee;
    position = this->mean_position + kalman_gains[j]*(*measurements[j] - myys[j]) + this->As[j]*second_part_shiftee;
  }
}
