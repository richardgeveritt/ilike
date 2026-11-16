#include "gaussian_smc_adaptor.h"
#include "utils.h"
#include "smc_output.h"
#include "ensemble_kalman_output.h"
#include "matrix_parameter_estimator.h"
#include "vector_parameter_estimator.h"

namespace ilike
{
GaussianSMCAdaptor::GaussianSMCAdaptor()
:SMCAdaptor()
{
}

GaussianSMCAdaptor::~GaussianSMCAdaptor()
{
  for (auto i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  
  for (auto i=this->covariance_estimators.begin();
       i!=this->covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
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
  
  for (auto i=this->mean_estimators.begin();
       i!=this->mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
  }
  this->mean_estimators.clear();
  
  for (auto i=this->covariance_estimators.begin();
       i!=this->covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      (*i).reset();
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
  for (auto i=another.mean_estimators.begin();
       i!=another.mean_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->mean_estimators.emplace_back((*i)->duplicate());
    else
      this->mean_estimators.emplace_back();
  }
  
  this->covariance_estimators.resize(0);
  this->covariance_estimators.reserve(another.covariance_estimators.size());
  for (auto i=another.covariance_estimators.begin();
       i!=another.covariance_estimators.end();
       ++i)
  {
    if (*i!=NULL)
      this->covariance_estimators.emplace_back((*i)->duplicate());
    else
      this->covariance_estimators.emplace_back();
  }
}

SMCAdaptor* GaussianSMCAdaptor::duplicate() const
{
  return( new GaussianSMCAdaptor(*this));
}

void GaussianSMCAdaptor::smc_adapt(SMCOutput* current_state)
{
  // reads the generation's particles one by one: not available for a
  // generation stored flat (dev/plan.md 2.6a; no adaptor is set from R)
  if (current_state->back().is_flat())
    Rcpp::stop("GaussianSMCAdaptor::smc_adapt: not available for a generation stored flat.");
  // the estimators take normalised log weights (they exponentiate them)
  arma::colvec wt = current_state->back().normalised_log_weights;
  // the estimators take the particles as plain pointers (not owned)
  std::vector<MoveOutput*> points;
  for (auto i = current_state->back().particles.begin(); i != current_state->back().particles.end(); ++i)
    points.push_back(i->get());
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    if (this->mean_estimators[i]!=NULL)
    {
      this->mean_estimators[i]->fit(variable_names[i],
                                    points,
                                    wt);
      this->gaussian_info_pointers[i]->get_mean() = this->mean_estimators[i]->estimated;
    }
    
    if (this->covariance_estimators[i]!=NULL)
    {
      this->covariance_estimators[i]->fit(variable_names[i],
                                          points,
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
  // the estimators take the ensemble members as plain pointers (not owned)
  std::vector<MoveOutput*> points;
  points.reserve(current_state->back().members.size());
  for (auto i = current_state->back().members.begin(); i != current_state->back().members.end(); ++i)
    points.push_back(i->get());
  for (size_t i=0; i<this->variable_names.size(); ++i)
  {
    if (this->mean_estimators[i]!=NULL)
    {
      this->mean_estimators[i]->fit(variable_names[i],
                                    points,
                                    wt);
      this->gaussian_info_pointers[i]->get_mean() = this->mean_estimators[i]->estimated;
    }
    
    if (this->covariance_estimators[i]!=NULL)
    {
      this->covariance_estimators[i]->fit(variable_names[i],
                                          points,
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
}
