#include <iterator>
#include "langevin_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"
#include "gradient_estimator.h"
#include "gradient_estimator_output.h"

LangevinProposalKernel::LangevinProposalKernel()
  :ProposalKernel()
{
  this->unused_variables_kept = true;
  this->index = NULL;
  this->gradient_estimator = NULL;
}

LangevinProposalKernel::~LangevinProposalKernel()
{
  if (this->gradient_estimator!=NULL)
    delete this->gradient_estimator;
  
  if (this->index!=NULL)
    delete index;
}

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in)
{
  this->unused_variables_kept = true;
  
  this->gradient_estimator = NULL;
  this->index = NULL;
  
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                                                           const std::vector<arma::mat> &covariances_in)
{
  this->unused_variables_kept = true;
  
  this->gradient_estimator = NULL;
  this->index = NULL;
  
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
}

LangevinProposalKernel::LangevinProposalKernel(const LangevinProposalKernel &another)
  :ProposalKernel(another)
{
  this->make_copy(another);
}

void LangevinProposalKernel::operator=(const LangevinProposalKernel &another)
{
  if(this == &another)
    return;
  
  if (this->gradient_estimator!=NULL)
    delete this->gradient_estimator;
  
  if (this->index!=NULL)
    delete index;

  ProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* LangevinProposalKernel::duplicate() const
{
  return( new LangevinProposalKernel(*this));
}

ProposalKernel* LangevinProposalKernel::proposal_kernel_duplicate() const
{
  return( new LangevinProposalKernel(*this));
}

void LangevinProposalKernel::make_copy(const LangevinProposalKernel &another)
{
  this->unused_variables_kept = another.unused_variables_kept;
  this->proposal_info = another.proposal_info;
  if (another.gradient_estimator!=NULL)
    this->gradient_estimator = another.gradient_estimator->duplicate();
  else
    this->gradient_estimator = NULL;
  
  if (another.index!=NULL)
    this->index = another.index->duplicate();
  else
    this->index = NULL;
}

double LangevinProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                        Particle &old_particle) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = old_particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_vector(i->first),
                              mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

double LangevinProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
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
    arma::colvec mean = old_particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_vector(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

double LangevinProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                  Particle &old_particle) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  /*
   double output = 0.0;
   for (size_t i=0; i<this->covariances.size(); ++i)
   {
   output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
   old_particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],old_particle,conditioned_on_parameters))/2.0,
   this->invs[i],
   this->logdets[i]);
   }
   return output;
   */
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = old_particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_vector(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

double LangevinProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                  Particle &old_particle,
                                                                  const Parameters &conditioned_on_parameters) const
{
  GradientEstimatorOutput* estimator = old_particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  /*
  double output = 0.0;
  for (size_t i=0; i<this->covariances.size(); ++i)
  {
    output = output + dmvnorm(proposed_particle.parameters.get_vector(this->variable_names[i]),
                              old_particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],old_particle,conditioned_on_parameters))/2.0,
                              this->invs[i],
                              this->logdets[i]);
  }
  return output;
   */
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = old_particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_vector(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

Parameters LangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                            Particle &particle) const
{
  /*
  Parameters output = particle.parameters;
  for (size_t i=0; i<this->covariances.size(); ++i)
  {
    output[this->variable_names[i]] = rmvnorm(rng,
                                              particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->get_gradient_of_log(this->variable_names[i],particle))/2.0,
                                              this->chols[i],
                                              true);
  }
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters LangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                            Particle &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  /*
  Parameters output = particle.parameters;
  for (size_t i=0; i<this->covariances.size(); ++i)
  {
    output[this->variable_names[i]] = rmvnorm(rng,
                                              particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->get_gradient_of_log(this->variable_names[i],particle,conditioned_on_parameters))/2.0,
                                              this->chols[i],
                                              true);
  }
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      Particle &particle) const
{
  /*
  Parameters output = particle.parameters;
  for (size_t i=0; i<this->covariances.size(); ++i)
  {
    output[this->variable_names[i]] = rmvnorm(rng,
                                              particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],particle))/2.0,
                                              this->chols[i],
                                              true);
  }
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,particle))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                            Particle &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  /*
  Parameters output = particle.parameters;
  for (size_t i=0; i<this->covariances.size(); ++i)
  {
    output[this->variable_names[i]] = rmvnorm(rng,
                                              particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],particle,conditioned_on_parameters))/2.0,
                                              this->chols[i],
                                              true);
  }
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_vector(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const std::string &variable,
                                                      Particle &particle) const
{
  /*
  Parameters output = particle.parameters;
  size_t i = [std::distance(this->variable_names.begin(),
                            std::find(this->variable_names.begin(),
                                      this->variable_names.end(),
                                      variable))];
  output[this->variable_names[i]] = rmvnorm(rng,
                                            particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],particle))/2.0,
                                            this->chols[i],
                                            true);
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  arma::colvec mean = particle.move_parameters->get_vector(found->first);
  double scale = found->second.get_double_scale();
  //double dim = double(mean.n_rows);
  output[found->first] = rmvnorm(rng,
                             mean + scale*found->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(found->first,this->index,particle))/2.0,
                             sqrt(scale)*found->second.get_chol(),
                             true);
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const std::string &variable,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  /*
  Parameters output = particle.parameters;
  size_t i = [std::distance(this->variable_names.begin(),
                            std::find(this->variable_names.begin(),
                                      this->variable_names.end(),
                                      variable))];
  output[this->variable_names[i]] = rmvnorm(rng,
                                            particle.parameters.get_vector(this->variable_names[i]) + this->covariances[i]*arma::vectorise(this->gradient_estimator->subsample_get_gradient_of_log(this->variable_names[i],particle,conditioned_on_parameters))/2.0,
                                            this->chols[i],
                                            true);
  return output;
  */
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  arma::colvec mean = particle.move_parameters->get_vector(found->first);
  double scale = found->second.get_double_scale();
  //double dim = double(mean.n_rows);
  output[found->first] = rmvnorm(rng,
                             mean + scale*found->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(found->first,this->index,particle,conditioned_on_parameters))/2.0,
                             sqrt(scale)*found->second.get_chol(),
                             true);
  return output;
}

arma::mat LangevinProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         Particle &proposed_particle,
                                                         Particle &old_particle)
{
  throw std::runtime_error("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat LangevinProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         Particle &proposed_particle,
                                                         Particle &old_particle,
                                                         const Parameters &conditioned_on_parameters)
{
  throw std::runtime_error("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat LangevinProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                   Particle &proposed_particle,
                                                                   Particle &old_particle,
                                                                   const Parameters &conditioned_on_parameters)
{
  throw std::runtime_error("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}
