#include "adjustment_ensemble_shifter.h"
#include "ensemble_factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble.h"

// follows https://arxiv.org/abs/2006.02941

AdjustmentEnsembleShifter::AdjustmentEnsembleShifter()
  :EnsembleShifter()
{
}

AdjustmentEnsembleShifter::~AdjustmentEnsembleShifter()
{
  
}

//Copy constructor for the AdjustmentEnsembleShifter class.
AdjustmentEnsembleShifter::AdjustmentEnsembleShifter(const AdjustmentEnsembleShifter &another)
  :EnsembleShifter(another)
{
  this->make_copy(another);
}

void AdjustmentEnsembleShifter::operator=(const AdjustmentEnsembleShifter &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  EnsembleShifter::operator=(another);
  this->make_copy(another);
}

EnsembleShifter* AdjustmentEnsembleShifter::duplicate() const
{
  return( new AdjustmentEnsembleShifter(*this));
}

void AdjustmentEnsembleShifter::make_copy(const AdjustmentEnsembleShifter &another)
{
  this->Zf = another.Zf;
  //this->Ginv = another.Ginv;
  //this->Ftranspose = another.Ftranspose;
  this->mean_position = another.mean_position;
  //this->Vs = another.Vs;
  this->As = another.As;
  
  this->P = another.P;
  this->Vtranspose = another.Vtranspose;
  this->Dhathalf = another.Dhathalf;
  this->Yhats = another.Yhats;
}

void AdjustmentEnsembleShifter::setup(Ensemble* ensemble,
                                      double inverse_incremental_temperature)
{
  this->mean_position = arma::mean(ensemble->packed_members,0).as_col();
  this->Zf = (1.0/sqrt(ensemble->size()-1.0))*(ensemble->packed_members.each_row() - this->mean_position.t()).t();
  
  arma::vec Dhathalf_vec;
  arma::svd(this->P,Dhathalf_vec,this->Vtranspose,this->Zf);
  //this->Ftranspose = this->Ftranspose.t();
  this->Dhathalf = arma::mat(this->Zf.n_rows,this->Zf.n_cols);
  
  /*
  for (size_t i=0; i<this->Zf.n_cols; ++i)
  {
    this->Dhathalf.col(i) =
  }
  */
  this->Dhathalf.diag() = Dhathalf_vec;
  
  this->Yhats.clear();
  this->Yhats.reserve(ensemble->Cxys.size());
  for (size_t i=0;
       i<ensemble->Cxys.size();
       ++i)
  {
    this->Yhats.push_back((1.0/sqrt(ensemble->size()-1.0))*(ensemble->packed_measurement_states[i].each_row()-arma::conv_to<arma::rowvec>::from(ensemble->myys[i])));
  }
  
  this->As = ensemble->ensemble_factors->get_adjustments(this->Zf,
                                                         this->Dhathalf,
                                                         this->P,
                                                         this->Vtranspose,
                                                         this->Yhats,
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

void AdjustmentEnsembleShifter::shift(const EnsembleFactorVariables* ensemble_factor_variables,
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
