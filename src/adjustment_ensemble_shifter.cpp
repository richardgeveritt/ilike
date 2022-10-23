#include "adjustment_ensemble_shifter.h"
#include "ensemble_factor_variables.h"
#include "ensemble.h"

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
  this->Ginv = another.Ginv;
  this->Ftranspose = another.Ftranspose;
  this->mean_position = another.mean_position;
  this->Vs = another.Vs;
}

void AdjustmentEnsembleShifter::setup(Ensemble* ensemble)
{
  this->mean_position = arma::mean(ensemble->packed_members,0);
  this->Zf = (1.0/(ensemble->size()-1.0))*(ensemble->packed_members.each_row() - this->mean_position).t();
  
  arma::mat U;
  arma::vec G;
  arma::svd_econ(this->Ftranspose,G,U,this->Zf);
  this->Ftranspose = this->Ftranspose.t();
  this->Ginv.diag() = G;
  this->Ginv = arma::pinv(this->Ginv);
  
  this->Vs.clear();
  this->Vs.reserve(ensemble->Cxys.size());
  for (size_t i=0;
       i<ensemble->Cxys.size();
       ++i)
  {
    arma::rowvec mean_meas_state = arma::mean(ensemble->packed_measurement_states[i],0);
    this->Vs.push_back(ensemble->packed_measurement_states[i].each_row()-mean_meas_state);
  }
}

void AdjustmentEnsembleShifter::shift(const EnsembleFactorVariables* ensemble_factor_variables,
                                      arma::colvec &position,
                                      const std::vector<arma::mat> &Cxys,
                                      const std::vector<arma::mat> &Cyys,
                                      double inverse_incremental_temperature) const
{
  std::vector<arma::colvec> shift_terms = ensemble_factor_variables->get_deterministic_shifts();
  
  std::vector<arma::mat> kalman_gains = ensemble_factor_variables->get_kalman_gains(Cxys,
                                                                                    Cyys,
                                                                                    inverse_incremental_temperature);
  
  std::vector<arma::colvec*> measurements = ensemble_factor_variables->get_measurements();
  
  std::vector<arma::mat> As = ensemble_factor_variables->get_adjustments(this->Zf,
                                                                         this->Ginv,
                                                                         this->Ftranspose,
                                                                         this->Vs,
                                                                         inverse_incremental_temperature);
  
  //std::vector<arma::colvec> first_part_of_shift;
  //first_part_of_shift.reserve(Cxys.size());
  
  arma::colvec second_part_shiftee = position - this->mean_position;
  position = this->mean_position;
  
  for (size_t j=0;
       j<Cxys.size();
       ++j)
  {
    position = position + kalman_gains[j]*(*measurements[j] - shift_terms[j]) + As[j]*second_part_shiftee;
  }
}
