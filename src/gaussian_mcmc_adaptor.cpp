#include "gaussian_mcmc_adaptor.h"
#include "utils.h"
#include "smc_output.h"
#include "gaussian_recursive_parameter_estimator.h"
#include "vector_recursive_parameter_estimator.h"
#include "scale_recursive_parameter_estimator.h"

namespace ilike
{
GaussianMCMCAdaptor::GaussianMCMCAdaptor()
:MCMCAdaptor()
{
}

GaussianMCMCAdaptor::~GaussianMCMCAdaptor()
{
  for (std::vector<VectorRecursiveParameterEstimator*>::iterator i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  for (std::vector<GaussianRecursiveParameterEstimator*>::iterator i=this->gaussian_estimators.begin();
       i!=this->gaussian_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  
  for (std::vector<ScaleRecursiveParameterEstimator*>::iterator i=this->scale_estimators.begin();
       i!=this->scale_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

GaussianMCMCAdaptor::GaussianMCMCAdaptor(const GaussianMCMCAdaptor &another)
:MCMCAdaptor(another)
{
  this->make_copy(another);
}

void GaussianMCMCAdaptor::operator=(const GaussianMCMCAdaptor &another)
{
  if(this == &another)
    return;
  
  for (std::vector<VectorRecursiveParameterEstimator*>::iterator i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->mean_estimators.clear();
  
  for (std::vector<GaussianRecursiveParameterEstimator*>::iterator i=this->gaussian_estimators.begin();
       i!=this->gaussian_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->gaussian_estimators.clear();
  
  for (std::vector<ScaleRecursiveParameterEstimator*>::iterator i=this->scale_estimators.begin();
       i!=this->scale_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
  this->scale_estimators.clear();
  
  MCMCAdaptor::operator=(another);
  this->make_copy(another);
}

void GaussianMCMCAdaptor::make_copy(const GaussianMCMCAdaptor &another)
{
  this->variable_names = another.variable_names;
  this->gaussian_info_pointers = another.gaussian_info_pointers;
  this->initial_proposal_info = another.initial_proposal_info;
  
  this->mean_estimators.resize(0);
  this->mean_estimators.reserve(another.mean_estimators.size());
  for (std::vector<VectorRecursiveParameterEstimator*>::const_iterator i=another.mean_estimators.begin();
       i!=another.mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->mean_estimators.push_back((*i)->vector_duplicate());
    else
      this->mean_estimators.push_back(NULL);
  }
  
  this->gaussian_estimators.resize(0);
  this->gaussian_estimators.reserve(another.gaussian_estimators.size());
  for (std::vector<GaussianRecursiveParameterEstimator*>::const_iterator i=another.gaussian_estimators.begin();
       i!=another.gaussian_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->gaussian_estimators.push_back((*i)->gaussian_duplicate());
    else
      this->gaussian_estimators.push_back(NULL);
  }
  
  this->scale_estimators.resize(0);
  this->scale_estimators.reserve(another.scale_estimators.size());
  for (std::vector<ScaleRecursiveParameterEstimator*>::const_iterator i=another.scale_estimators.begin();
       i!=another.scale_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->scale_estimators.push_back((*i)->scale_duplicate());
    else
      this->scale_estimators.push_back(NULL);
  }
}

MCMCAdaptor* GaussianMCMCAdaptor::duplicate() const
{
  return( new GaussianMCMCAdaptor(*this));
}

void GaussianMCMCAdaptor::specific_mcmc_adapt(const Particle &latest_particle,
                                              size_t iteration_counter)
{
  Rcpp::stop("GaussianMCMCAdaptor::specific_mcmc_adapt - need to project particle into the correct space each time. This needs to be done in ther update step.");
  
  // need to know what proposal we are adapting for, since need to know if we are in transformed space
  // need to set move parameters in base adapt?
  // actually, do this here
  
  // also for this gaussian adaptor, might want to have some initial GaussianPropopalInfo and use this in the early iterations
  
  // particle stores acceptance info for purposes of adapting, gets copied when we propose a new particle, gets wiped clean when we adapt: do this at end of composite adapt
  
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    if (this->mean_estimators[i]!=NULL)
    {
      this->mean_estimators[i]->update(variable_names[i],
                                       latest_particle,
                                       iteration_counter,
                                       this->proposal);
      this->gaussian_info_pointers[i]->get_mean() = this->mean_estimators[i]->estimated;
    }
    
    if (this->gaussian_estimators[i]!=NULL)
    {
      this->gaussian_estimators[i]->update(variable_names[i],
                                           latest_particle,
                                           iteration_counter,
                                           this->proposal);
      *this->gaussian_info_pointers[i] = this->gaussian_estimators[i]->estimated;
      this->gaussian_info_pointers[i]->set_covariance_info();
    }
    
    if (this->scale_estimators[i]!=NULL)
    {
      this->scale_estimators[i]->update(variable_names[i],
                                        latest_particle,
                                        iteration_counter,
                                        this->proposal);
      this->gaussian_info_pointers[i]->get_scale() = this->scale_estimators[i]->estimated;
    }
  }
  
  
  // Find sample covariance of current particle set.
  //arma::mat matrix_particles = current_state->back().get_most_recent_matrix_particles();
  //arma::colvec wt = exp(current_state->back().normalised_log_weights);
  
  //arma::mat sample_covariance = cov_wt(matrix_particles,wt);
  //size_t dimension = sample_covariance.n_rows;
  
  // needs sorting
  // Set specified variable in parameter_to_adapt to have scaled covariance.
  //parameter_to_adapt[this->variance_name] = this->scale(dimension)*sample_covariance;
}
}
