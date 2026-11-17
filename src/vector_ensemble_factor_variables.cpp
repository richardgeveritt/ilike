#include "vector_ensemble_factor_variables.h"
#include "vector_ensemble_factors.h"
#include "measurement_covariance_estimator_output.h"
#include "ensemble_factors.h"
#include "index.h"

namespace ilike
{
VectorEnsembleFactorVariables::VectorEnsembleFactorVariables()
:EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs.resize(0);
  this->vector_ensemble_factors = NULL;
}

VectorEnsembleFactorVariables::VectorEnsembleFactorVariables(const VectorEnsembleFactors* vector_ensemble_factors_in,
                                                             const std::vector<MeasurementCovarianceEstimatorOutput*> &measurement_covariance_estimators_in)
:EnsembleFactorVariables()
{
  this->measurement_covariance_estimator_outputs.reserve(measurement_covariance_estimators_in.size());
  for (auto k : measurement_covariance_estimators_in)
    this->measurement_covariance_estimator_outputs.emplace_back(k);
  this->vector_ensemble_factors = vector_ensemble_factors_in;
}

VectorEnsembleFactorVariables::~VectorEnsembleFactorVariables()
{
  for (auto i=this->measurement_covariance_estimator_outputs.begin();
       i!=this->measurement_covariance_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
}

//Copy constructor for the VectorEnsembleFactorVariables class.
VectorEnsembleFactorVariables::VectorEnsembleFactorVariables(const VectorEnsembleFactorVariables &another)
:EnsembleFactorVariables(another)
{
  this->make_copy(another);
}

void VectorEnsembleFactorVariables::operator=(const VectorEnsembleFactorVariables &another)
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

EnsembleFactorVariables* VectorEnsembleFactorVariables::duplicate() const
{
  return( new VectorEnsembleFactorVariables(*this));
}

void VectorEnsembleFactorVariables::make_copy(const VectorEnsembleFactorVariables &another)
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
  this->vector_ensemble_factors = another.vector_ensemble_factors;
}

std::vector<arma::rowvec> VectorEnsembleFactorVariables::get_measurement_states_for_covariance() const
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

std::vector<arma::colvec> VectorEnsembleFactorVariables::get_shifts(double inverse_incremental_temperature) const
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

std::vector<arma::colvec> VectorEnsembleFactorVariables::get_deterministic_shifts() const
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

double VectorEnsembleFactorVariables::evaluate_ensemble_likelihood_ratios(const Index* index,
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

double VectorEnsembleFactorVariables::subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
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

double VectorEnsembleFactorVariables::evaluate_likelihoods(const Index* index) const
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

double VectorEnsembleFactorVariables::subsample_evaluate_likelihoods(const Index* index) const
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

const EnsembleFactors* VectorEnsembleFactorVariables::get_ensemble_factors() const
{
  return this->vector_ensemble_factors;
}

void VectorEnsembleFactorVariables::forget_you_were_already_written_to_file()
{
}

void VectorEnsembleFactorVariables::write_to_file(const std::string &directory_name,
                                                  const std::string &index) const
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    std::string factor_directory_name = directory_name + "/ensemble_factor" + toString(i);
    this->measurement_covariance_estimator_outputs[i]->write_to_file(factor_directory_name,index);
  }
}

void VectorEnsembleFactorVariables::close_ofstreams()
{
  for (size_t i=0;
       i<this->measurement_covariance_estimator_outputs.size();
       ++i)
  {
    this->measurement_covariance_estimator_outputs[i]->close_ofstreams();
  }
}

}
