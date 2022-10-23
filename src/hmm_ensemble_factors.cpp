#include "hmm_ensemble_factors.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "hmm_ensemble_factor_variables.h"
#include "proposal_kernel.h"

HMMEnsembleFactors::HMMEnsembleFactors()
  :EnsembleFactors()
{
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimator_temp_data.resize(0);
  this->transition_kernel = NULL;
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
  
  for (std::vector<Data*>::iterator i=this->measurement_covariance_estimator_temp_data.begin();
       i!=this->measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
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
  
  for (std::vector<Data*>::iterator i=this->measurement_covariance_estimator_temp_data.begin();
       i!=this->measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->measurement_covariance_estimator_temp_data.clear();
  
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
  
  this->measurement_covariance_estimator_temp_data.resize(0);
  for (std::vector<Data*>::const_iterator i=another.measurement_covariance_estimator_temp_data.begin();
       i!=another.measurement_covariance_estimator_temp_data.end();
       ++i)
  {
    if (*i!=NULL)
      this->measurement_covariance_estimator_temp_data.push_back((*i)->duplicate());
    else
      this->measurement_covariance_estimator_temp_data.push_back(NULL);
  }
}

void HMMEnsembleFactors::set_data(const Index* index)
{
  if (this->measurement_covariance_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
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
  for (auto i=this->measurement_covariance_estimators.begin();
       i!=this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->change_data(subsetted_data);
  }
}

EnsembleFactorVariables* HMMEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters)
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

arma::colvec HMMEnsembleFactors::get_measurements()
{
  return measurement_covariance_estimator_temp_data[0]->get_vector(this->measurement_names);
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
                                       const std::vector<arma::mat> &Cyys)
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

/*
std::vector<arma::mat> VectorEnsembleFactors::get_kalman_gains(const std::vector<arma::mat> &Cxys,
                                                               const std::vector<arma::mat> &Cyys,
                                                               double inverse_incremental_temperature)
{
  std::vector<arma::mat> kalman_gains;
  kalman_gains.reserve(this->measurement_covariance_estimators.size);
  for (size_t i=0;
       i<this->measurement_covariance_estimators.size();
       ++i)
  {
    this->measurement_covariance_estimators[i]->get_kalman_gain(Cxys[i],
                                                                Cyys[i],
                                                                inverse_incremental_temperature);
  }
}
*/
