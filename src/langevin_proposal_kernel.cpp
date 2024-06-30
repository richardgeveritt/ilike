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
  //this->unused_variables_kept = true;
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

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  //this->unused_variables_kept = true;
  
  this->gradient_estimator = gradient_estimator_in;
  this->gradient_estimator->set_proposal(this);
  this->index = NULL;
  
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

LangevinProposalKernel::LangevinProposalKernel(const std::string &variable_name_in,
                                               const arma::mat &covariance_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  //this->unused_variables_kept = true;
  
  this->gradient_estimator = gradient_estimator_in;
  this->gradient_estimator->set_proposal(this);
  this->index = NULL;
  
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
}

LangevinProposalKernel::LangevinProposalKernel(const std::string &variable_name_in,
                                               const double &sd_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  //this->unused_variables_kept = true;
  
  this->gradient_estimator = gradient_estimator_in;
  this->gradient_estimator->set_proposal(this);
  this->index = NULL;
  
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
}

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                                               const std::vector<arma::mat> &covariances_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  //this->unused_variables_kept = true;
  
  this->gradient_estimator = gradient_estimator_in;
  this->gradient_estimator->set_proposal(this);
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
  //this->unused_variables_kept = another.unused_variables_kept;
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

double LangevinProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                        const Particle &old_particle) const
{
  GradientEstimatorOutput* gradient_estimator_output = old_particle.get_gradient_estimator_output(this);
  
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = old_particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(i->first),
                              mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index,old_particle))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
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
    arma::colvec mean = old_particle.move_parameters->get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_colvec(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}
*/

double LangevinProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                  const Particle &old_particle) const
{
  GradientEstimatorOutput* gradient_estimator_output = old_particle.get_gradient_estimator_output(this);
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
    arma::colvec mean = old_particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index,old_particle))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
double LangevinProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
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
    arma::colvec mean = old_particle.move_parameters->get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.move_parameters->get_colvec(i->first),
                                            mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,old_particle,conditioned_on_parameters))/2.0,
                                            (1.0/sqrt(i->second.get_double_scale()))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}
*/

Parameters LangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                            const Particle &particle) const
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
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = particle.get_transformed_parameters(this);
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm_using_chol(rng,
                                          mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index,particle))/2.0,
                                          sqrt(scale)*i->second.get_chol());
  }
  return output;
}

/*
Parameters LangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                            Particle &particle,
                                            const Parameters &conditioned_on_parameters) const
{
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                     this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_colvec(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}
*/

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const Particle &particle) const
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
  
  Rcpp::stop("LangevinProposalKernel::subsample_simulate - not written yet.");
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = particle.get_transformed_parameters(this);
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm_using_chol(rng,
                                          mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(i->first,this->index,particle))/2.0,
                                          sqrt(scale)*i->second.get_chol());
  }
  return output;
}

/*
Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.move_parameters->get_colvec(i->first);
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm(rng,
                               mean + scale*i->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(i->first,this->index,particle,conditioned_on_parameters))/2.0,
                               sqrt(scale)*i->second.get_chol(),
                               true);
  }
  return output;
}
*/

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const std::string &variable,
                                                      const Particle &particle) const
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
  
  Rcpp::stop("LangevinProposalKernel::subsample_simulate - not written yet.");
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = particle.get_transformed_parameters(this);
  arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(found->first);
  double scale = found->second.get_double_scale();
  //double dim = double(mean.n_rows);
  output[found->first] = rmvnorm_using_chol(rng,
                                            mean + scale*found->second.get_covariance()*arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(found->first,this->index,particle))/2.0,
                                            sqrt(scale)*found->second.get_chol());
  return output;
}

/*
Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const std::string &variable,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  
  GradientEstimatorOutput* estimator = particle.initialise_gradient_estimator_output(this,
                                                                                         this->gradient_estimator);
  
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  if (this->unused_variables_kept)
    output = *particle.move_parameters;
  arma::colvec mean = particle.move_parameters->get_colvec(found->first);
  double scale = found->second.get_double_scale();
  //double dim = double(mean.n_rows);
  output[found->first] = rmvnorm(rng,
                             mean + scale*found->second.get_covariance()*arma::vectorise(estimator->subsample_get_gradient_of_log(found->first,this->index,particle,conditioned_on_parameters))/2.0,
                             sqrt(scale)*found->second.get_chol(),
                             true);
  return output;
}
*/

arma::mat LangevinProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           const Particle &proposed_particle,
                                                           const Particle &old_particle)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat LangevinProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                         Particle &proposed_particle,
                                                         Particle &old_particle,
                                                         const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat LangevinProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                     const Particle &proposed_particle,
                                                                     const Particle &old_particle)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat LangevinProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                   Particle &proposed_particle,
                                                                   Particle &old_particle,
                                                                   const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

void LangevinProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

GradientEstimatorOutput* LangevinProposalKernel::simulate_gradient_estimator_output() const
{
  GradientEstimatorOutput* current_output = gradient_estimator->initialise();
  current_output->simulate_auxiliary_variables();
  return current_output;
}

std::vector<const ProposalKernel*> LangevinProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void LangevinProposalKernel::set_index(Index* index_in)
{
  this->index = index_in;
}

void LangevinProposalKernel::set_index_if_null(Index* index_in)
{
  if (this->index==NULL)
    this->index = index_in;
}

bool LangevinProposalKernel::can_be_evaluated() const
{
  return true;
}

void LangevinProposalKernel::set_data(Data* data_in)
{
  
}
