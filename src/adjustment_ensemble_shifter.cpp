#include "adjustment_ensemble_shifter.h"
#include "ensemble_factors.h"
#include "ensemble_factor_variables.h"
#include "ensemble.h"

namespace ilike
{
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
  this->mean_position = another.mean_position;
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
  
  // the compact SVD Zf = P Dhat^(1/2) V^T, with only the r non-zero singular
  // values (Anderson 2001; Huang et al., arXiv-2204.04386, the EAKI): P is
  // d x r, Dhat^(1/2) r x r, and the matrix stored is V^T (r x J), as the
  // estimators' get_adjustment use it (V^T M V, an r x r matrix). (Before:
  // the full SVD, and V stored as if it were V^T, so the adjustment used
  // V M V^T over all J directions and collapsed the ensemble.)
  arma::mat U_full;
  arma::vec singular_values;
  arma::mat V_full;
  arma::svd_econ(U_full, singular_values, V_full, this->Zf);
  double tolerance = (singular_values.n_elem > 0) ?
    singular_values.max() * double(std::max(this->Zf.n_rows, this->Zf.n_cols)) * arma::datum::eps : 0.0;
  arma::uword rank = arma::accu(singular_values > tolerance);
  if (rank == 0)
    Rcpp::stop("AdjustmentEnsembleShifter: the ensemble has no spread.");
  this->P = U_full.cols(0, rank - 1);
  this->Vtranspose = V_full.cols(0, rank - 1).t();
  this->Dhathalf = arma::diagmat(singular_values.subvec(0, rank - 1));
  
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
  
}

void AdjustmentEnsembleShifter::shift(const EnsembleFactorVariables* ensemble_factor_variables,
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
