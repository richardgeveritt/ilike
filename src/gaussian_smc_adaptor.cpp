#include "gaussian_smc_adaptor.h"
#include "utils.h"
#include "smc_output.h"
#include "ensemble_kalman_output.h"
#include "matrix_parameter_estimator.h"
#include "vector_parameter_estimator.h"

GaussianSMCAdaptor::GaussianSMCAdaptor()
  :SMCAdaptor()
{
}

GaussianSMCAdaptor::~GaussianSMCAdaptor()
{
  for (std::vector<VectorParameterEstimator*>::iterator i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  for (std::vector<MatrixParameterEstimator*>::iterator i=this->covariance_estimators.begin();
       i!=this->covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

GaussianSMCAdaptor::GaussianSMCAdaptor(const GaussianSMCAdaptor &another)
  :SMCAdaptor(another)
{
  this->make_copy(another);
}

void GaussianSMCAdaptor::operator=(const GaussianSMCAdaptor &another)
{
  if(this == &another)
    return;
  
  for (std::vector<VectorParameterEstimator*>::iterator i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->mean_estimators.clear();
  
  for (std::vector<MatrixParameterEstimator*>::iterator i=this->covariance_estimators.begin();
       i!=this->covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->covariance_estimators.clear();

  SMCAdaptor::operator=(another);
  this->make_copy(another);
}

void GaussianSMCAdaptor::make_copy(const GaussianSMCAdaptor &another)
{
  this->variable_names = another.variable_names;
  this->scales = another.scales;
  this->gaussian_info_pointers = another.gaussian_info_pointers;
  
  this->mean_estimators.resize(0);
  this->mean_estimators.reserve(another.mean_estimators.size());
  for (std::vector<VectorParameterEstimator*>::const_iterator i=another.mean_estimators.begin();
       i!=another.mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->mean_estimators.push_back((*i)->duplicate());
    else
      this->mean_estimators.push_back(NULL);
  }
  
  this->covariance_estimators.resize(0);
  this->covariance_estimators.reserve(another.covariance_estimators.size());
  for (std::vector<MatrixParameterEstimator*>::const_iterator i=another.covariance_estimators.begin();
       i!=another.covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->covariance_estimators.push_back((*i)->duplicate());
    else
      this->covariance_estimators.push_back(NULL);
  }
}

SMCAdaptor* GaussianSMCAdaptor::duplicate() const
{
  return( new GaussianSMCAdaptor(*this));
}

void GaussianSMCAdaptor::smc_adapt(SMCOutput* current_state)
{
  arma::colvec wt = exp(current_state->back().normalised_log_weights);
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    if (this->mean_estimators[i]!=NULL)
    {
      this->mean_estimators[i]->fit(variable_names[i],
                                    current_state->back().particles,
                                    wt);
      this->gaussian_info_pointers[i]->get_mean() = this->mean_estimators[i]->estimated;
    }
    
    if (this->covariance_estimators[i]!=NULL)
    {
      this->covariance_estimators[i]->fit(variable_names[i],
                                          current_state->back().particles,
                                          wt);
      this->gaussian_info_pointers[i]->get_covariance() = covariance_estimators[i]->estimated;
    }
  }
  
  // Find sample covariance of current particle set.
  //arma::mat matrix_particles = current_state->back().get_most_recent_matrix_particles();
  
  //arma::mat sample_covariance = cov_wt(matrix_particles,wt);
  
  // needs sorting
  // Set specified variable in parameter_to_adapt to have scaled covariance.
  //parameter_to_adapt[this->variance_name] = this->scale(dimension)*sample_covariance;
}

void GaussianSMCAdaptor::ensemble_adapt(EnsembleKalmanOutput* current_state)
{
  arma::colvec wt(current_state->back().size());
  wt.fill(-log(double(current_state->back().size())));
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    if (this->mean_estimators[i]!=NULL)
    {
      this->mean_estimators[i]->fit(variable_names[i],
                                    current_state->back().members,
                                    wt);
      this->gaussian_info_pointers[i]->get_mean() = this->mean_estimators[i]->estimated;
    }
    
    if (this->covariance_estimators[i]!=NULL)
    {
      this->covariance_estimators[i]->fit(variable_names[i],
                                          current_state->back().members,
                                          wt);
      this->gaussian_info_pointers[i]->get_covariance() = covariance_estimators[i]->estimated;
    }
  }
  
  // Find sample covariance of current particle set.
  //arma::mat matrix_particles = current_state->back().get_most_recent_matrix_particles();
  
  //arma::mat sample_covariance = cov_wt(matrix_particles,wt);
  
  // needs sorting
  // Set specified variable in parameter_to_adapt to have scaled covariance.
  //parameter_to_adapt[this->variance_name] = this->scale(dimension)*sample_covariance;
}
