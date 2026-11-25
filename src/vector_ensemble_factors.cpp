#include "vector_ensemble_factors.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "vector_factor_variables.h"
#include "vector_ensemble_factor_variables.h"
#include "parameters.h"
#include "index.h"
#include "ensemble.h"

namespace ilike
{
VectorEnsembleFactors::VectorEnsembleFactors()
:EnsembleFactors()
{
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimator_temp_data.resize(0);
}

VectorEnsembleFactors::VectorEnsembleFactors(const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in)
:EnsembleFactors()
{
  this->measurement_covariance_estimators.reserve(measurement_covariance_estimators_in.size());
  for (auto k : measurement_covariance_estimators_in)
    this->measurement_covariance_estimators.emplace_back(k);
}

VectorEnsembleFactors::~VectorEnsembleFactors()
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  
}

//Copy constructor for the VectorEnsembleFactors class.
VectorEnsembleFactors::VectorEnsembleFactors(const VectorEnsembleFactors &another)
:EnsembleFactors(another)
{
  this->make_copy(another);
}

void VectorEnsembleFactors::operator=(const VectorEnsembleFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  this->measurement_covariance_estimators.clear();

  EnsembleFactors::operator=(another);
  this->make_copy(another);
}

EnsembleFactors* VectorEnsembleFactors::duplicate() const
{
  return( new VectorEnsembleFactors(*this));
}

void VectorEnsembleFactors::make_copy(const VectorEnsembleFactors &another)
{
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimators.reserve(another.measurement_covariance_estimators.size());
  for (auto i=another.measurement_covariance_estimators.begin();
       i!=another.measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimators.emplace_back((*i)->duplicate());
    else
      this->measurement_covariance_estimators.emplace_back();
  }
  
  this->measurement_covariance_estimator_temp_data = another.measurement_covariance_estimator_temp_data;
}

void VectorEnsembleFactors::set_data(const Index* index)
{
  if (this->measurement_covariance_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  
  // assumption that all llhd_estimators point to the same data
  Data* all_data = this->measurement_covariance_estimators[0]->get_data();
  
  if (this->measurement_covariance_estimator_temp_data.size()==0)
  {
    this->measurement_covariance_estimator_temp_data.push_back(std::make_shared<Data>(all_data->rows(indices)));
  }
  else
  {
    this->measurement_covariance_estimator_temp_data[0] = std::make_shared<Data>(all_data->rows(indices));
  }
  
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->change_data(this->measurement_covariance_estimator_temp_data[0]);
  }
}

std::vector<arma::colvec*> VectorEnsembleFactors::get_measurements()
{
  std::vector<arma::colvec*> measurements;
  measurements.reserve(this->measurement_covariance_estimators.size());
  
  for (auto i=this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement_pointer());
  }
  
  return measurements;
}

EnsembleFactorVariables* VectorEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  for (auto i=this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->simulate(simulated_parameters);
  }
  
  return new VectorEnsembleFactorVariables(this,
                                           outputs);
}

EnsembleFactorVariables* VectorEnsembleFactors::subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  for (auto i=this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->subsample_simulate(simulated_parameters);
  }
  
  return new VectorEnsembleFactorVariables(this,
                                           outputs);
}

bool VectorEnsembleFactors::need_Cxx() const
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if ((*i)->need_Cxx()==true)
    {
      return true;
      break;
    }
  }
  return false;
}

void VectorEnsembleFactors::find_Cygivenx(const arma::mat &inv_Cxx,
                                          const std::vector<arma::mat> &Cxys,
                                          const std::vector<arma::mat> &Cyys) const
{
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    this->measurement_covariance_estimators[i]->find_Cygivenx(inv_Cxx,
                                                              Cxys[i],
                                                              Cyys[i]);
  }
}

std::vector<arma::mat> VectorEnsembleFactors::get_adjustments(const arma::mat &Zf,
                                                              const arma::mat &Dhathalf,
                                                              const arma::mat &P,
                                                              const arma::mat &Vtranspose,
                                                              const std::vector<arma::mat> &Yhats,
                                                              double inverse_incremental_temperature) const
{
  std::vector<arma::mat> adjustments;
  adjustments.reserve(this->measurement_covariance_estimators.size());
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    adjustments.push_back(this->measurement_covariance_estimators[i]->get_adjustment(Zf,
                                                                                     Dhathalf,
                                                                                     P,
                                                                                     Vtranspose,
                                                                                     Yhats[i],
                                                                                     inverse_incremental_temperature));
  }
  return adjustments;
}

std::vector<arma::mat> VectorEnsembleFactors::get_sqrt_adjustments(const std::vector<arma::mat> &Cxys,
                                                                   const std::vector<arma::mat> &Cyys,
                                                                   double inverse_incremental_temperature) const
{
  std::vector<arma::mat> adjustments;
  adjustments.reserve(this->measurement_covariance_estimators.size());
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    adjustments.push_back(this->measurement_covariance_estimators[i]->get_sqrt_adjustment(Cxys[i],
                                                                                          Cyys[i],
                                                                                          inverse_incremental_temperature));
  }
  return adjustments;
}

double VectorEnsembleFactors::get_incremental_likelihood(Ensemble* ensemble) const
{
  double inverse_incremental_temperature = 1.0;///(this->temperature - this->previous_temperature);
  
  double llhd = 0.0;
  ensemble->kalman_gains.clear();
  ensemble->kalman_gains.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],
                                                                                                                                          inverse_incremental_temperature);
    
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    llhd = llhd + dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance);
  }
  
  return llhd;
}

double VectorEnsembleFactors::get_mc_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                   double inverse_incremental_temperature) const
{
  
  double llhd = 0.0;
  ensemble->kalman_gains.clear();
  ensemble->kalman_gains.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],
                                                                                                                                          inverse_incremental_temperature);
    
    arma::mat Cygivenx = this->measurement_covariance_estimators[i]->get_Cygivenx();
    double d = Cygivenx.n_rows;
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    
    // _without_sympd version used to due possible numerical isses with Sigma matrices for bad parameters
    llhd = llhd + (d/2.0)*log(inverse_incremental_temperature) + (d/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1/2.0)*(1.0-(1.0/inverse_incremental_temperature))*arma::log_det_sympd(Cygivenx) + dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance);
  }
  
  return llhd;
}

double VectorEnsembleFactors::get_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                   double inverse_incremental_temperature) const
{
  
  double llhd = 0.0;
  ensemble->kalman_gains.clear();
  ensemble->kalman_gains.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],
                                                                                                                                          inverse_incremental_temperature);
    
    arma::mat Cygivenx = this->measurement_covariance_estimators[i]->get_Cygivenx();
    double d = Cygivenx.n_rows;
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());

    double c = (d/2.0)*log(inverse_incremental_temperature) + (d/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1.0/2.0)*(1.0-(1.0/inverse_incremental_temperature))*arma::log_det_sympd(Cygivenx);

    // _without_sympd version used to due possible numerical isses with Sigma matrices for bad parameters
    llhd = llhd + c + dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance);
  }
  
  return llhd;
}

double VectorEnsembleFactors::get_unbiased_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                            double inverse_incremental_temperature) const
{
  double llhd = 0.0;
  ensemble->kalman_gains.clear();
  ensemble->kalman_gains.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],
                                                                                                                                          inverse_incremental_temperature);
    
    arma::mat Cygivenx = this->measurement_covariance_estimators[i]->get_Cygivenx();
    double d = Cygivenx.n_rows;
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    
    size_t n = ensemble->size();
    if (n<d+4)
    {
      Rcpp::stop("VectorEnsembleFactors::get_unbiased_inversion_incremental_likelihood - to use the unbiased estimator ensemble size needs to be greater than d+3.");
    }
    
    llhd = llhd + (d/2.0)*log(inverse_incremental_temperature) + (d/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1/2.0)*(1.0-(1.0/inverse_incremental_temperature))*arma::log_det_sympd(Cygivenx) +  dmvnorm_estimated_params(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance,n);
  }
  
  return llhd;
}

void VectorEnsembleFactors:: get_path1_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                        std::vector<double> &log_measurement_likelihood_means,
                                                                        double temperature,
                                                                        double multiplier) const
{
  log_measurement_likelihood_means.clear();
  log_measurement_likelihood_means.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat Cygivenx = this->measurement_covariance_estimators[i]->get_Cygivenx();
    
    size_t n = ensemble->size();
    
    arma::colvec log_likelihoods(n);
    for (size_t j=0; j<n; ++j)
    {
      log_likelihoods[j] = dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(), ensemble->partially_packed_members_col[j], Cygivenx);
    }

    log_measurement_likelihood_means.push_back(multiplier*arma::mean(log_likelihoods));

  }
  
}

void VectorEnsembleFactors::get_path2_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                       std::vector<double> &log_measurement_likelihood_means,
                                                                       std::vector<double> &log_measurement_likelihood_variances) const
{
  
  log_measurement_likelihood_means.clear();
  log_measurement_likelihood_means.reserve(this->measurement_covariance_estimators.size());
  
  log_measurement_likelihood_variances.clear();
  log_measurement_likelihood_variances.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    
    arma::mat Cygivenx = this->measurement_covariance_estimators[i]->get_Cygivenx();
    
    size_t n = ensemble->size();
    arma::colvec log_likelihoods(n);
    for (size_t j=0; j<n; ++j)
    {
      log_likelihoods[j] = dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(), ensemble->partially_packed_members_col[j], Cygivenx);
    }
    log_measurement_likelihood_means.push_back(arma::mean(log_likelihoods));
    log_measurement_likelihood_variances.push_back(arma::var(log_likelihoods));
  }
}

void VectorEnsembleFactors::calculate_kalman_gains(Ensemble* ensemble,
                                                   double inverse_incremental_temperature) const
{
  ensemble->kalman_gains.clear();
  ensemble->kalman_gains.reserve(this->measurement_covariance_estimators.size());
  
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    arma::mat unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],
                                                                                                                                          inverse_incremental_temperature);
    
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
  }
  
}

void VectorEnsembleFactors::setup()
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void VectorEnsembleFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}

void VectorEnsembleFactors::precompute_gaussian_covariance(double inverse_incremental_temperature,
                                                           std::vector<arma::mat> &inv_sigma_precomps,
                                                           std::vector<double> &log_det_precomps) const
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    arma::mat inv_sigma_precomp;
    double log_det_precomp;
    (*i)->precompute_gaussian_covariance(inverse_incremental_temperature,
                                         inv_sigma_precomp,
                                         log_det_precomp);
    inv_sigma_precomps.push_back(inv_sigma_precomp);
    log_det_precomps.push_back(log_det_precomp);
  }
}


// (C280) member i's weight from its rows, as its outputs' ratios add up
double VectorEnsembleFactors::compact_ensemble_likelihood_ratios(const Index* index,
                                               const Ensemble* ensemble,
                                               size_t i,
                                               double inverse_incremental_temperature,
                                               const std::vector<arma::mat> &inv_sigma_precomps,
                                               const std::vector<double> &log_det_precomps) const
{
  double result = 0.0;
  for (auto j=index->begin(); j!=index->end(); ++j)
    result = result + this->measurement_covariance_estimators[*j]->ensemble_likelihood_ratio(arma::colvec(ensemble->compact_measurement_states[*j].row(i).t()),
                                                                                            inverse_incremental_temperature,
                                                                                            inv_sigma_precomps[*j],
                                                                                            log_det_precomps[*j]);
  return result;
}

bool VectorEnsembleFactors::supports_compact_members() const
{
  if (this->measurement_covariance_estimators.empty())
    return false;
  for (auto i = this->measurement_covariance_estimators.begin(); i != this->measurement_covariance_estimators.end(); ++i)
    if ((*i == NULL) || !(*i)->supports_compact_members())
      return false;
  return true;
}

arma::colvec VectorEnsembleFactors::compact_measurement_state(size_t j, const Parameters &parameters)
{
  return this->measurement_covariance_estimators[j]->compact_measurement_state(parameters);
}

arma::colvec VectorEnsembleFactors::compact_noise(size_t j)
{
  return this->measurement_covariance_estimators[j]->compact_noise();
}

void VectorEnsembleFactors::compact_simulate(size_t j, const Parameters &parameters, arma::colvec &measurement_state, arma::colvec &noise)
{
  this->measurement_covariance_estimators[j]->compact_simulate(parameters, measurement_state, noise);
}

arma::colvec VectorEnsembleFactors::compact_shift(size_t j, const arma::colvec &measurement_state, const arma::colvec &noise, double inverse_incremental_temperature) const
{
  return this->measurement_covariance_estimators[j]->compact_shift(measurement_state, noise, inverse_incremental_temperature);
}
}
