#include <math.h>
#include <iterator>
#include "barker_dynamics_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"
#include "gradient_estimator.h"
#include "independent_proposal_kernel.h"
#include "gradient_estimator_output.h"

BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel()
  :ProposalKernel()
{
  this->gradient_estimator = NULL;
  this->proposal_simulate = NULL;
  this->index = NULL;
}

BarkerDynamicsProposalKernel::~BarkerDynamicsProposalKernel()
{
  if (this->gradient_estimator!=NULL)
    delete this->gradient_estimator;
  
  if (this->proposal_simulate!=NULL)
    delete this->proposal_simulate;
  
  if (this->index!=NULL)
    delete index;
}

BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in)
{
  this->gradient_estimator = NULL;
  this->proposal_simulate = NULL;
  this->index = NULL;
  
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const std::vector<std::string> &variable_names_in,
                                                           const std::vector<arma::mat> &covariances_in)
{
  this->gradient_estimator = NULL;
  this->index = NULL;
  
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
}

BarkerDynamicsProposalKernel::BarkerDynamicsProposalKernel(const BarkerDynamicsProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void BarkerDynamicsProposalKernel::operator=(const BarkerDynamicsProposalKernel &another)
{
  if(this == &another)
    return;
  
  if (this->gradient_estimator!=NULL)
    delete this->gradient_estimator;
  
  if (this->proposal_simulate!=NULL)
    delete this->proposal_simulate;
  
  if (this->index!=NULL)
    delete index;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* BarkerDynamicsProposalKernel::duplicate() const
{
  return( new BarkerDynamicsProposalKernel(*this));
}

ProposalKernel* BarkerDynamicsProposalKernel::proposal_kernel_duplicate() const
{
  return( new BarkerDynamicsProposalKernel(*this));
}

void BarkerDynamicsProposalKernel::make_copy(const BarkerDynamicsProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  if (another.gradient_estimator!=NULL)
    this->gradient_estimator = another.gradient_estimator->duplicate();
  else
    this->gradient_estimator = NULL;
  if (another.proposal_simulate!=NULL)
    this->proposal_simulate = another.proposal_simulate->independent_proposal_kernel_duplicate();
  else
    this->proposal_simulate = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;

}

double BarkerDynamicsProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                              Particle &old_particle) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                        this->gradient_estimator);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle));
    
    arma::mat proposed_mat = proposed_particle.parameters[i->first];
    arma::mat old_mat = old_particle.parameters[i->first];
    arma::colvec current_proposed = arma::vectorise(proposed_mat);
    arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();
    
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
    }
  }
  return output;
}

double BarkerDynamicsProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                              Particle &old_particle,
                                                              const Parameters &conditioned_on_parameters) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters));
    
    arma::mat proposed_mat = proposed_particle.parameters[i->first];
    arma::mat old_mat = old_particle.parameters[i->first];
    arma::colvec current_proposed = arma::vectorise(proposed_mat);
    arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();
    
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
    }
  }
  return output;
}

double BarkerDynamicsProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                        Particle &old_particle) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle));
    
    arma::mat proposed_mat = proposed_particle.parameters[i->first];
    arma::mat old_mat = old_particle.parameters[i->first];
    arma::colvec current_proposed = arma::vectorise(proposed_mat);
    arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();
    
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
    }
  }
  return output;
}

double BarkerDynamicsProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                        Particle &old_particle,
                                                                        const Parameters &conditioned_on_parameters) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle));
    
    arma::mat proposed_mat = proposed_particle.parameters[i->first];
    arma::mat old_mat = old_particle.parameters[i->first];
    arma::colvec current_proposed = arma::vectorise(proposed_mat);
    arma::rowvec z = (current_proposed - arma::vectorise(old_mat))*i->second.get_inv_chol();
    
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      output = output - log1p(exp(-z[j]*current_chol_gradient_prod[j]));
    }
  }
  return output;
}

Parameters BarkerDynamicsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                  Particle &particle) const
{
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  Parameters proposed = this->proposal_simulate->simulate(rng,particle);
  
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle));
    
    arma::mat initial_proposed = proposed[i->first];
    arma::rowvec current_proposed = arma::vectorise(initial_proposed);
    
    arma::rowvec log_unifs = log(multiple_runif(rng,current_proposed.size()));
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
      if (log_unifs[j]>log_prob)
      {
        current_proposed[j] = -current_proposed[j];
      }
    }
    proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * sqrt(i->second.get_double_scale())*i->second.get_chol(),initial_proposed.n_rows,initial_proposed.n_cols);
  }
  return proposed;
}

Parameters BarkerDynamicsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                  Particle &particle,
                                                  const Parameters &conditioned_on_parameters) const
{
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  Parameters proposed = this->proposal_simulate->simulate(rng,
                                                          particle,
                                                          conditioned_on_parameters);
  
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = sqrt(i->second.get_double_scale())*i->second.get_chol() * arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters));
    
    arma::mat initial_proposed = proposed[i->first];
    arma::rowvec current_proposed = arma::vectorise(initial_proposed);
    
    arma::rowvec log_unifs = log(multiple_runif(rng,current_proposed.size()));
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
      if (log_unifs[j]>log_prob)
      {
        current_proposed[j] = -current_proposed[j];
      }
    }
    proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * sqrt(i->second.get_double_scale())*i->second.get_chol(),initial_proposed.n_rows,initial_proposed.n_cols);
  }
  return proposed;
}

Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            Particle &particle) const
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  /*
  Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
                                                particle);
  
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = i->second.chol * arma::vectorise(this->subsample_gradient_estimator->get_gradient_of_log(i->first,particle));
    
    arma::mat initial_proposed = proposed[i->first];
    arma::rowvec current_proposed = arma::vectorise(initial_proposed);
    
    arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
      if (log_unifs[j]>log_prob)
      {
        current_proposed[j] = -current_proposed[j];
      }
    }
    proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * i->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
  }
  return proposed;
  */
}
                                   
Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            Particle &particle,
                                                            const Parameters &conditioned_on_parameters) const
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  /*
  Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
                                                particle,
                                                conditioned_on_parameters);
  
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_chol_gradient_prod = i->second.chol * arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(i->first,particle,conditioned_on_parameters));
    
    arma::mat initial_proposed = proposed[i->first];
    arma::rowvec current_proposed = arma::vectorise(initial_proposed);
    
    arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
    for (size_t j=0; j<current_proposed.size(); ++j)
    {
      double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
      if (log_unifs[j]>log_prob)
      {
        current_proposed[j] = -current_proposed[j];
      }
    }
    proposed[i->first] = particle.parameters[i->first] + arma::reshape(current_proposed * i->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
  }
  return proposed;
  */
}

Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            const std::string &variable,
                                                            Particle &particle) const
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  /*
  Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
                                                          variable,
                                                          particle);
  
  auto found = this->proposal_info.find(variable);
  
  arma::colvec current_chol_gradient_prod = found->second.chol * arma::vectorise(this->subsample_gradient_estimator->get_gradient_of_log(found->first,particle));
  
  arma::mat initial_proposed = proposed[found->first];
  arma::rowvec current_proposed = arma::vectorise(initial_proposed);
  
  arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
  for (size_t j=0; j<current_proposed.size(); ++j)
  {
    double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
    if (log_unifs[j]>log_prob)
    {
      current_proposed[j] = -current_proposed[j];
    }
  }
  proposed[found->first] = particle.parameters[found->first] + arma::reshape(current_proposed * found->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
  return proposed;
  */
}

Parameters BarkerDynamicsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                            const std::string &variable,
                                                       Particle &particle,
                                                       const Parameters &conditioned_on_parameters) const
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::subsample_simulate - not written yet.");
  /*
  Parameters proposed = this->proposal_simulate->subsample_simulate(rng,
                                                          variable,
                                                          particle,
                                                          conditioned_on_parameters);
  
  auto found = this->proposal_info.find(variable);
  
  arma::colvec current_chol_gradient_prod = found->second.chol * arma::vectorise(this->subsample_gradient_estimator->get_gradient_of_log(found->first,particle,conditioned_on_parameters));
  
  arma::mat initial_proposed = proposed[found->first];
  arma::rowvec current_proposed = arma::vectorise(initial_proposed);
  
  arma::rowvec log_unifs = log(runif(rng,current_proposed.size()));
  for (size_t j=0; j<current_proposed.size(); ++j)
  {
    double log_prob = -log1p(exp(-current_proposed[j]*current_chol_gradient_prod[j]));
    if (log_unifs[j]>log_prob)
    {
      current_proposed[j] = -current_proposed[j];
    }
  }
  proposed[found->first] = particle.parameters[found->first] + arma::reshape(current_proposed * found->second.chol,initial_proposed.n_rows,initial_proposed.n_cols);
  return proposed;
  */
}

arma::mat BarkerDynamicsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                   Particle &proposed_particle,
                                   Particle &old_particle)
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat BarkerDynamicsProposalKernel::specific_gradient_of_log(const std::string &variable,
                                   Particle &proposed_particle,
                                   Particle &old_particle,
                                   const Parameters &conditioned_on_parameters)
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
}

//virtual arma::mat specific_subsample_gradient_of_log(const std::string &variable,
//                                                     Particle &proposed_particle,
//                                                     Particle &old_particle)=0;
arma::mat BarkerDynamicsProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                             Particle &proposed_particle,
                                             Particle &old_particle,
                                             const Parameters &conditioned_on_parameters)
{
  throw std::runtime_error("BarkerDynamicsProposalKernel::specific_gradient_of_log - not written yet.");
}
