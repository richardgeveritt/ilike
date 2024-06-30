#include "hmm_ensemble_factors.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "hmm_ensemble_factor_variables.h"
#include "proposal_kernel.h"
#include "ensemble.h"

HMMEnsembleFactors::HMMEnsembleFactors()
  :EnsembleFactors()
{
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimator_temp_data.resize(0);
  this->transition_kernel = NULL;
}

HMMEnsembleFactors::HMMEnsembleFactors(ProposalKernel* transition_kernel_in,
                                       const std::vector<MeasurementCovarianceEstimator*> &measurement_covariance_estimators_in)
:EnsembleFactors()
{
  this->transition_kernel = transition_kernel_in;
  this->measurement_covariance_estimators = measurement_covariance_estimators_in;
  this->measurement_covariance_estimator_temp_data.resize(0);
}

HMMEnsembleFactors::~HMMEnsembleFactors()
{
  for (std::vector<MeasurementCovarianceEstimator*>::iterator i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  /*
  for (auto i=this->measurement_covariance_estimator_temp_data.begin();
       i!=this->measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  */
  
  if (this->transition_kernel!=NULL)
    delete this->transition_kernel;
}

//Copy constructor for the HMMEnsembleFactors class.
HMMEnsembleFactors::HMMEnsembleFactors(const HMMEnsembleFactors &another)
  :EnsembleFactors(another)
{
  this->make_copy(another);
}

void HMMEnsembleFactors::operator=(const HMMEnsembleFactors &another)
{
  if(this == &another){ //if a==a
    return;
  }
  
  if (this->transition_kernel!=NULL)
    delete this->transition_kernel;

  for (std::vector<MeasurementCovarianceEstimator*>::iterator i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->measurement_covariance_estimators.clear();
  
  /*
  for (auto i=this->measurement_covariance_estimator_temp_data.begin();
       i!=this->measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->measurement_covariance_estimator_temp_data.clear();
  */
  
  EnsembleFactors::operator=(another);
  this->make_copy(another);
}

EnsembleFactors* HMMEnsembleFactors::duplicate() const
{
  return( new HMMEnsembleFactors(*this));
}

void HMMEnsembleFactors::make_copy(const HMMEnsembleFactors &another)
{
  if (another.transition_kernel!=NULL)
    this->transition_kernel = another.transition_kernel->proposal_kernel_duplicate();
  else
    this->transition_kernel = NULL;
  
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimators.reserve(another.measurement_covariance_estimators.size());
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i=another.measurement_covariance_estimators.begin();
       i!=another.measurement_covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimators.push_back((*i)->duplicate());
    else
      this->measurement_covariance_estimators.push_back(NULL);
  }
  
  this->measurement_covariance_estimator_temp_data = another.measurement_covariance_estimator_temp_data;
  /*
  this->measurement_covariance_estimator_temp_data.resize(0);
  for (auto i=another.measurement_covariance_estimator_temp_data.begin();
       i!=another.measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimator_temp_data.push_back((*i)->duplicate());
    else
      this->measurement_covariance_estimator_temp_data.push_back(NULL);
  }
  */
}

void HMMEnsembleFactors::set_data(const Index* index)
{
  if (this->measurement_covariance_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  
  /*
  for (int i=0; i<indices.size(); ++i)
  {
  
  }
  */
  
  /*
  Data* subsetted_data = new Data();
  
  // assumption that all llhd_estimators point to the same data
  Data* all_data = this->measurement_covariance_estimators[0]->get_data();
  
  for (auto i=all_data->vector_begin();
       i!=all_data->vector_end();
       ++i)
  {
    (*subsetted_data)[i->first] = (*all_data)[i->first].rows(indices);
  }
   
  if (this->measurement_covariance_estimator_temp_data[0]!=NULL)
    delete this->measurement_covariance_estimator_temp_data[0];
   
  this->measurement_covariance_estimator_temp_data[0] = subsetted_data;
  */
  
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

EnsembleFactorVariables* HMMEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->simulate(simulated_parameters);
  }
  
  return new HMMEnsembleFactorVariables(this,
                                        outputs);
}

/*
EnsembleFactorVariables* HMMEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
                                                                                const Parameters &conditioned_on_parameters)
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  Parameters all_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(all_parameters));
    outputs.back()->simulate(all_parameters);
  }
  
  return new HMMEnsembleFactorVariables(this,
                                        outputs);
}
*/

EnsembleFactorVariables* HMMEnsembleFactors::subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters) const
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(simulated_parameters));
    outputs.back()->subsample_simulate(simulated_parameters);
  }
  
  return new HMMEnsembleFactorVariables(this,
                                        outputs);
}

/*
EnsembleFactorVariables* HMMEnsembleFactors::subsample_simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
                                                                                const Parameters &conditioned_on_parameters)
{
  std::vector<MeasurementCovarianceEstimatorOutput*> outputs;
  outputs.reserve(this->measurement_covariance_estimators.size());
  
  Parameters all_parameters = simulated_parameters.merge(conditioned_on_parameters);
  
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    outputs.push_back((*i)->initialise(all_parameters));
    outputs.back()->subsample_simulate(all_parameters);
  }
  
  return new HMMEnsembleFactorVariables(this,
                                        outputs);
}
*/

std::vector<arma::colvec*> HMMEnsembleFactors::get_measurements()
{
  std::vector<arma::colvec*> measurements;
  measurements.reserve(this->measurement_covariance_estimators.size());
  
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    measurements.push_back((*i)->get_measurement_pointer());
  }
  
  return measurements;
  //return measurement_covariance_estimator_temp_data[0]->get_colvec(this->measurement_names);
}

/*
std::vector<arma::mat> HMMEnsembleFactors::get_measurement_covariances()
{
  std::vector<arma::mat> output;
  output.reserve(this->measurement_covariance_estimators.size());
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    output.push_back((*i)->get_measurement_covariamce());
  }
  return output;
}
 
std::vector<arma::mat> HMMEnsembleFactors::get_measurement_covariances(const Parameters &conditioned_on_parameters)
{
  std::vector<arma::mat> output;
  output.reserve(this->measurement_covariance_estimators.size());
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    output.push_back((*i)->get_measurement_covariamce(conditioned_on_parameters));
  }
  return output;
}
*/

bool HMMEnsembleFactors::need_Cxx() const
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

void HMMEnsembleFactors::find_Cygivenx(const arma::mat &inv_Cxx,
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

std::vector<arma::mat> HMMEnsembleFactors::get_adjustments(const arma::mat &Zf,
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

std::vector<arma::mat> HMMEnsembleFactors::get_sqrt_adjustments(const std::vector<arma::mat> &Cxys,
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

double HMMEnsembleFactors::get_incremental_likelihood(Ensemble* ensemble) const
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
    //arma::mat non_tempered_unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],1.0);
    
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    llhd = llhd + dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance);
  }
  
  return llhd;
}

double HMMEnsembleFactors::get_inversion_incremental_likelihood(Ensemble* ensemble,
                                                                double inverse_incremental_temperature) const
{
  //double inverse_incremental_temperature = 1.0/(this->temperature - this->previous_temperature);
  
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
    //arma::mat non_tempered_unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],1.0);
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    llhd = llhd + (d/2.0)*log(inverse_incremental_temperature) + (d/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1/2.0)*(1.0-(1.0/inverse_incremental_temperature))*arma::log_det_sympd(Cygivenx) + dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance);
  }
  
  return llhd;
}

double HMMEnsembleFactors::get_unbiased_inversion_incremental_likelihood(Ensemble* ensemble,
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
    //arma::mat non_tempered_unconditional_measurement_covariance = this->measurement_covariance_estimators[i]->get_unconditional_measurement_covariance(ensemble->Cyys[i],1.0);
    ensemble->kalman_gains.push_back(ensemble->Cxys[i]*unconditional_measurement_covariance.i());
    
    size_t n = ensemble->size();
    if (n<d+4)
    {
      Rcpp::stop("HMMEnsembleFactors::get_unbiased_inversion_incremental_likelihood - to use the unbiased estimator ensemble size needs to be greater than d+3.");
    }
    
    llhd = llhd + (d/2.0)*log(inverse_incremental_temperature) + (d/2.0)*(1.0-(1.0/inverse_incremental_temperature))*log(2.0*M_PI) + (1/2.0)*(1.0-(1.0/inverse_incremental_temperature))*arma::log_det_sympd(Cygivenx) +  dmvnorm_estimated_params(*this->measurement_covariance_estimators[i]->get_measurement_pointer(),ensemble->myys[i],unconditional_measurement_covariance,n);
  }
  
  return llhd;
}

void HMMEnsembleFactors::get_path1_inversion_incremental_likelihood(Ensemble* ensemble,
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

    //double mean_for_this_term = 0.0;
    size_t n = ensemble->size();
    arma::colvec log_likelihoods(n);
    for (size_t j=0; j<n; ++j)
    {
      log_likelihoods[j] = dmvnorm(*this->measurement_covariance_estimators[i]->get_measurement_pointer(), ensemble->partially_packed_members_col[j], Cygivenx);
    }
    log_measurement_likelihood_means.push_back(arma::mean(log_likelihoods));
  }
}

void HMMEnsembleFactors::get_path2_inversion_incremental_likelihood(Ensemble* ensemble,
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

void HMMEnsembleFactors::calculate_kalman_gains(Ensemble* ensemble,
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

void HMMEnsembleFactors::setup()
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup();
  }
}

void HMMEnsembleFactors::setup(const Parameters &conditioned_on_parameters)
{
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->setup(conditioned_on_parameters);
  }
}

void HMMEnsembleFactors::precompute_gaussian_covariance(double inverse_incremental_temperature,
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
