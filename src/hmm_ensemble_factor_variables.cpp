#include "hmm_ensemble_factor_variables.h"
#include "hmm_ensemble_factors.h"
#include "index.h"
#include "measurement_covariance_estimator_output.h"

namespace ilike
{
HMMEnsembleFactorVariables::HMMEnsembleFactorVariables()
:EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs.resize(0);
  this->hmm_ensemble_factors = NULL;
}

HMMEnsembleFactorVariables::~HMMEnsembleFactorVariables()
{
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
}

HMMEnsembleFactorVariables::HMMEnsembleFactorVariables(const HMMEnsembleFactors* hmm_ensemble_factors_in,
                                                       const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimator_outputs_in)
: EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs.reserve(measurement_covariance_estimator_outputs_in.size());
  for (auto k : measurement_covariance_estimator_outputs_in)
    this->measurement_covariance_estimator_outputs.emplace_back(k);
  this->hmm_ensemble_factors = hmm_ensemble_factors_in;
}

//Copy constructor for the HMMEnsembleFactorVariables class.
HMMEnsembleFactorVariables::HMMEnsembleFactorVariables(const HMMEnsembleFactorVariables &another)
:EnsembleFactorVariables(another)
{
  this->make_copy(another);
}

void HMMEnsembleFactorVariables::operator=(const HMMEnsembleFactorVariables &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  this->measurement_covariance_estimator_outputs.clear();
  
  EnsembleFactorVariables::operator=(another);
  this->make_copy(another);
}

EnsembleFactorVariables* HMMEnsembleFactorVariables::duplicate() const
{
  return( new HMMEnsembleFactorVariables(*this));
}

void HMMEnsembleFactorVariables::make_copy(const HMMEnsembleFactorVariables &another)
{
  this->measurement_covariance_estimator_outputs.resize(0);
  this->measurement_covariance_estimator_outputs.reserve(another.measurement_covariance_estimator_outputs.size());
  for (auto i=another.measurement_covariance_estimator_outputs.begin();
       i!=another.measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimator_outputs.emplace_back((*i)->duplicate());
    else
      this->measurement_covariance_estimator_outputs.emplace_back();
  }
  this->hmm_ensemble_factors = another.hmm_ensemble_factors;
}

std::vector<arma::rowvec> HMMEnsembleFactorVariables::get_measurement_states_for_covariance() const
{
  std::vector<arma::rowvec> measurements;
  measurements.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement_state_for_covariance());
  }
  return measurements;
}

std::vector<arma::colvec> HMMEnsembleFactorVariables::get_shifts(double inverse_incremental_temperature) const
{
  std::vector<arma::colvec> shifts;
  shifts.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    shifts.push_back((*i)->get_shift(inverse_incremental_temperature));
  }
  return shifts;
}

std::vector<arma::colvec> HMMEnsembleFactorVariables::get_deterministic_shifts() const
{
  std::vector<arma::colvec> shifts;
  shifts.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    shifts.push_back((*i)->get_deterministic_shift());
  }
  return shifts;
}

std::vector<arma::colvec*> HMMEnsembleFactorVariables::get_measurements() const
{
  std::vector<arma::colvec*> measurements;
  measurements.reserve(this->measurement_covariance_estimator_outputs.size());
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement());
  }
  return measurements;
}

double HMMEnsembleFactorVariables::evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                       double inverse_incremental_temperature,
                                                                       const std::vector<arma::mat> &inv_sigma_precomps,
                                                                       const std::vector<double> &log_det_precomps) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[*i]->evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature,
                                                                                                             inv_sigma_precomps[*i],
                                                                                                             log_det_precomps[*i]);
  }
  return result;
}

double HMMEnsembleFactorVariables::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                                                 double inverse_incremental_temperature,
                                                                                 const std::vector<arma::mat> &inv_sigma_precomps,
                                                                                 const std::vector<double> &log_det_precomps) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate_ensemble_likelihood_ratio(inverse_incremental_temperature,
                                                                                                                       inv_sigma_precomps[*i],
                                                                                                                       log_det_precomps[*i]);
  }
  return result;
}

double HMMEnsembleFactorVariables::evaluate_likelihoods(const Index* index) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[*i]->evaluate_likelihood();
  }
  return result;
}

double HMMEnsembleFactorVariables::subsample_evaluate_likelihoods(const Index* index) const
{
  double result = 0.0;
  for (auto i=index->begin();
       i!=index->end();
       ++i)
  {
    result = result + this->measurement_covariance_estimator_outputs[*i]->subsample_evaluate_likelihood();
  }
  return result;
}

const EnsembleFactors* HMMEnsembleFactorVariables::get_ensemble_factors() const
{
  return this->hmm_ensemble_factors;
}

void HMMEnsembleFactorVariables::forget_you_were_already_written_to_file()
{
}

void HMMEnsembleFactorVariables::write_to_file(const std::string &directory_name,
                                               const std::string &index) const
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    std::string factor_directory_name = directory_name + "/ensemble_factor" + toString(i);
    this->measurement_covariance_estimator_outputs[i]->write_to_file(factor_directory_name,
                                                                     index);
  }
}

void HMMEnsembleFactorVariables::close_ofstreams()
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    this->measurement_covariance_estimator_outputs[i]->close_ofstreams();
  }
}

// (C265) point at a copied filter's factors, and each measurement output at
// the estimator in the same place among its factors
bool HMMEnsembleFactorVariables::repoint_ensemble_factors(const EnsembleFactors* ensemble_factors_in)
{
  const HMMEnsembleFactors* factors = dynamic_cast<const HMMEnsembleFactors*>(ensemble_factors_in);
  if ((factors == NULL) || (factors->measurement_covariance_estimators.size() != this->measurement_covariance_estimator_outputs.size()))
    return false;
  this->hmm_ensemble_factors = factors;
  for (size_t i = 0; i < this->measurement_covariance_estimator_outputs.size(); ++i)
    if ((this->measurement_covariance_estimator_outputs[i] != NULL) &&
        !this->measurement_covariance_estimator_outputs[i]->repoint_estimator(factors->measurement_covariance_estimators[i].get()))
      return false;
  return true;
}

}
