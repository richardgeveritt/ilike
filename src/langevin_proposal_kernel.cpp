#include <iterator>
#include "langevin_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"
#include "gradient_estimator.h"
#include "gradient_estimator_output.h"

namespace ilike
{

LangevinProposalKernel::LangevinProposalKernel()
:ProposalKernel()
{
  this->index.reset();
  this->gradient_estimator.reset();
}

LangevinProposalKernel::~LangevinProposalKernel()
{
  this->gradient_estimator.reset();
  
  this->index.reset();
}

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  
  this->gradient_estimator.reset(gradient_estimator_in);
  this->gradient_estimator->set_proposal(this);
  this->index.reset();
  
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
  
  this->gradient_estimator.reset(gradient_estimator_in);
  this->gradient_estimator->set_proposal(this);
  this->index.reset();
  
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
}

LangevinProposalKernel::LangevinProposalKernel(const std::string &variable_name_in,
                                               const arma::mat &covariance_in,
                                               double scale_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  
  this->gradient_estimator.reset(gradient_estimator_in);
  this->gradient_estimator->set_proposal(this);
  this->index.reset();
  
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in,
                                                               scale_in);
}

LangevinProposalKernel::LangevinProposalKernel(const std::string &variable_name_in,
                                               const double &sd_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  
  this->gradient_estimator.reset(gradient_estimator_in);
  this->gradient_estimator->set_proposal(this);
  this->index.reset();
  
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
}

LangevinProposalKernel::LangevinProposalKernel(const std::vector<std::string> &variable_names_in,
                                               const std::vector<arma::mat> &covariances_in,
                                               GradientEstimator* gradient_estimator_in)
:ProposalKernel()
{
  
  this->gradient_estimator.reset(gradient_estimator_in);
  this->gradient_estimator->set_proposal(this);
  this->index.reset();
  
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
  
  this->gradient_estimator.reset();
  
  this->index.reset();
  
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
  this->proposal_info = another.proposal_info;
  if (another.gradient_estimator!=NULL)
    this->gradient_estimator.reset(another.gradient_estimator->duplicate());
  else
    this->gradient_estimator.reset();
  
  if (another.index!=NULL)
    this->index.reset(another.index->duplicate());
  else
    this->index.reset();
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
                                            mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index.get(),old_particle))/2.0,
                                            (1.0/i->second.get_double_scale())*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

double LangevinProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
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
                                            mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index.get(),old_particle))/2.0,
                                            (1.0/i->second.get_double_scale())*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

// simulate, written into the target (same random numbers and values)
void LangevinProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                           const Particle &particle,
                                           ProposalTarget &target) const
{
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    target.matrix(i->first) = rmvnorm_using_lower_chol(rng,
                                                       mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index.get(),particle))/2.0,
                                                       i->second.get_chol_lower(), sqrt(scale));
  }
}

Parameters LangevinProposalKernel::simulate(RandomNumberGenerator &rng,
                                            const Particle &particle) const
{
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    output[i->first] = rmvnorm_using_lower_chol(rng,
                                          mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->get_gradient_of_log(i->first,this->index.get(),particle))/2.0,
                                          i->second.get_chol_lower(), sqrt(scale));
  }
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const Particle &particle) const
{
  
  Rcpp::stop("LangevinProposalKernel::subsample_simulate - not written yet.");
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    output[i->first] = rmvnorm_using_lower_chol(rng,
                                          mean + scale*i->second.get_covariance()*arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(i->first,this->index.get(),particle))/2.0,
                                          i->second.get_chol_lower(), sqrt(scale));
  }
  return output;
}

Parameters LangevinProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      const std::string &variable,
                                                      const Particle &particle) const
{
  
  Rcpp::stop("LangevinProposalKernel::subsample_simulate - not written yet.");
  
  GradientEstimatorOutput* gradient_estimator_output = particle.get_gradient_estimator_output(this);
  
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  arma::colvec mean = particle.get_transformed_parameters(this).get_colvec(found->first);
  double scale = found->second.get_double_scale();
  output[found->first] = rmvnorm_using_lower_chol(rng,
                                            mean + scale*found->second.get_covariance()*arma::vectorise(gradient_estimator_output->subsample_get_gradient_of_log(found->first,this->index.get(),particle))/2.0,
                                            found->second.get_chol_lower(), sqrt(scale));
  return output;
}

arma::mat LangevinProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                           const Particle &proposed_particle,
                                                           const Particle &old_particle)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat LangevinProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                     const Particle &proposed_particle,
                                                                     const Particle &old_particle)
{
  Rcpp::stop("LangevinProposalKernel::specific_gradient_of_log - not written yet.");
}

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

// The kernel keeps its own private copy of the index: the pointer passed in is
// owned (and deleted) by the calling MCMC object, so storing it directly would
// lead to it being deleted twice.
void LangevinProposalKernel::set_index(Index* index_in)
{
  this->index.reset();
  if (index_in!=NULL)
    this->index.reset(index_in->duplicate());
  else
    this->index.reset();
}

void LangevinProposalKernel::set_index_if_null(Index* index_in)
{
  if ( (this->index==NULL) && (index_in!=NULL) )
    this->index.reset(index_in->duplicate());
}

bool LangevinProposalKernel::can_be_evaluated() const
{
  return true;
}

void LangevinProposalKernel::set_data(Data* data_in)
{
  
}

GaussianProposalInfo* LangevinProposalKernel::get_gaussian_proposal_info(const std::string &variable)
{
  auto found = this->proposal_info.find(variable);
  if (found==this->proposal_info.end())
    return NULL;
  return &found->second;
}

}
