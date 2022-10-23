#include "vector_ensemble_factors.h"
#include "measurement_covariance_estimator.h"
#include "measurement_covariance_estimator_output.h"
#include "vector_factor_variables.h"
#include "vector_ensemble_factor_variables.h"
#include "data.h"
#include "index.h"

VectorEnsembleFactors::VectorEnsembleFactors()
  :EnsembleFactors()
{
  this->measurement_covariance_estimators.resize(0);
  this->measurement_covariance_estimator_temp_data.resize(0);
}

VectorEnsembleFactors::~VectorEnsembleFactors()
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

EnsembleFactors* VectorEnsembleFactors::duplicate() const
{
  return( new VectorEnsembleFactors(*this));
}

void VectorEnsembleFactors::make_copy(const VectorEnsembleFactors &another)
{
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

void VectorEnsembleFactors::set_data(const Index* index)
{
  if (this->measurement_covariance_estimators.size()==0)
    return;
  
  arma::uvec indices = index->get_uvec();
  Data* subsetted_data = new Data();
  
  // assumption that all cov_estimators point to the same data
  Data* all_data = this->measurement_covariance_estimators[0]->get_data();
  
  for (auto i=this->measurement_names.begin();
       i!=this->measurement_names.end();
       ++i)
  {
    (*subsetted_data)[*i] = (*all_data)[*i].rows(indices);
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

arma::colvec VectorEnsembleFactors::get_measurements()
{
  return measurement_covariance_estimator_temp_data[0]->get_vector(this->measurement_names);
}

EnsembleFactorVariables* VectorEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters)
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
  
  return new VectorEnsembleFactorVariables(this,
                                           outputs);
}

EnsembleFactorVariables* VectorEnsembleFactors::simulate_ensemble_factor_variables(const Parameters &simulated_parameters,
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
    
    // call set params first
    // should be a proposalkernel? check for consistency
    
    //(*i)->set_parameters(conditioned_on_parameters);
    
    outputs.back()->simulate(all_parameters);
    // need to separate current state from params in this case
  }
  
  return new VectorEnsembleFactorVariables(this,
                                           outputs);
}

/*
std::vector<arma::mat> VectorEnsembleFactors::get_measurement_covariances()
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

std::vector<arma::mat> VectorEnsembleFactors::get_measurement_covariances(const Parameters &conditioned_on_parameters)
{
  std::vector<arma::mat> output;
  output.reserve(this->measurement_covariance_estimators.size());
  for (std::vector<MeasurementCovarianceEstimator*>::const_iterator i = this->measurement_covariance_estimators.begin();
       i != this->measurement_covariance_estimators.end();
       ++i)
  {
    (*i)->set_parameters(conditioned_on_parameters);
    output.push_back((*i)->get_measurement_covariamce());
  }
  return output;
}
*/

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
void VectorEnsembleFactors::find_measurement_covariances(EnsembleKalmanOutput* simulation)
{
  
}
*/
