#include <iterator>
#include "gaussian_random_walk_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel()
  :SymmetricProposalKernel()
{
  //this->unused_variables_kept = true;
}

GaussianRandomWalkProposalKernel::~GaussianRandomWalkProposalKernel()
{
}

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in)
  :SymmetricProposalKernel()
{
  //this->unused_variables_kept = true;
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = GaussianProposalInfo();
  }
}

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                   const std::vector<arma::mat> &covariances_in)
  :SymmetricProposalKernel()
{
  //this->unused_variables_kept = true;
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = GaussianProposalInfo(covariances_in[i]);
  }
}

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                   const double &sd_in)
:SymmetricProposalKernel()
{
  //this->unused_variables_kept = true;
  this->proposal_info[variable_name_in] = GaussianProposalInfo(sd_in);
}

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                   const arma::mat &covariance_in)
:SymmetricProposalKernel()
{
  //this->unused_variables_kept = true;
  this->proposal_info[variable_name_in] = GaussianProposalInfo(covariance_in);
}

GaussianRandomWalkProposalKernel::GaussianRandomWalkProposalKernel(const GaussianRandomWalkProposalKernel &another)
  :SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void GaussianRandomWalkProposalKernel::operator=(const GaussianRandomWalkProposalKernel &another)
{
  if(this == &another)
    return;

  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* GaussianRandomWalkProposalKernel::duplicate() const
{
  return( new GaussianRandomWalkProposalKernel(*this));
}

ProposalKernel* GaussianRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return( new GaussianRandomWalkProposalKernel(*this));
}

SymmetricProposalKernel* GaussianRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new GaussianRandomWalkProposalKernel(*this));
}

void GaussianRandomWalkProposalKernel::make_copy(const GaussianRandomWalkProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
  //this->unused_variables_kept = another.unused_variables_kept;
}

double GaussianRandomWalkProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                  const Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec mean = old_particle.get_transformed_parameters(this).get_colvec(i->first);
    double scale = i->second.get_double_scale();
    double dim = double(mean.n_rows);
    output = output + dmvnorm_using_precomp(proposed_particle.get_transformed_parameters(this).get_colvec(i->first),
                                            mean,
                                            (1.0/sqrt(scale))*i->second.get_inv(),
                                            dim*log(scale)+i->second.get_logdet());
  }
  return output;
}

/*
double GaussianRandomWalkProposalKernel::specific_evaluate_kernel(Particle &proposed_particle,
                                                                  Particle &old_particle,
                                                                  const Parameters &conditioned_on_parameters) const
{
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}
*/

double GaussianRandomWalkProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                            const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

/*
double GaussianRandomWalkProposalKernel::specific_subsample_evaluate_kernel(Particle &proposed_particle,
                                                                            Particle &old_particle,
                                                                            const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}
*/

void GaussianRandomWalkProposalKernel::set_covariance(const std::string &variable,
                                                      const arma::mat &covariance_in)
{
  this->proposal_info[variable].set_covariance(covariance_in);
}

arma::mat GaussianRandomWalkProposalKernel::get_inverse_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_inv();
}

arma::mat GaussianRandomWalkProposalKernel::get_covariance(const std::string &variable)
{
  return this->proposal_info[variable].get_covariance();
}

Parameters GaussianRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                      const Particle &particle) const
{
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    //arma::colvec mean = i->second.get_mean();
    double scale = i->second.get_double_scale();
    //double dim = double(mean.n_rows);
    output[i->first] = rmvnorm_using_chol(rng,
                                          particle.get_transformed_parameters(this).get_colvec(i->first),
                                          sqrt(scale)*i->second.get_chol());
  }
  return output;
}

/*
Parameters GaussianRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng, particle);
}
*/

Parameters GaussianRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}

/*
Parameters GaussianRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                      Particle &particle,
                                                      const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng, particle);
}
*/

Parameters GaussianRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const std::string &variable,
                                                                const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found = this->proposal_info.find(variable);
  
  Parameters output;
  //if (this->unused_variables_kept)
  //  output = *particle.move_parameters;
  output[variable] = rmvnorm_using_chol(rng,
                                        particle.get_transformed_parameters(this).get_colvec(variable),
                                        sqrt(found->second.get_double_scale())*found->second.get_chol());
  return output;
}

/*
Parameters GaussianRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                const std::string &variable,
                                                                Particle &particle,
                                                                const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->subsample_simulate(rng, variable, particle);
}
*/

arma::mat GaussianRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                     const Particle &proposed_particle,
                                                                     const Particle &old_particle)
{
  Rcpp::stop("GaussianRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat GaussianRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                    Particle &proposed_particle,
                                                                    Particle &old_particle,
                                                                    const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("GaussianRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

arma::mat GaussianRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                               const Particle &proposed_particle,
                                                                               const Particle &old_particle)
{
  Rcpp::stop("GaussianRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

/*
arma::mat GaussianRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                              Particle &proposed_particle,
                                                                              Particle &old_particle,
                                                                              const Parameters &conditioned_on_parameters)
{
  Rcpp::stop("GaussianRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}
*/

void GaussianRandomWalkProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

std::vector<std::string> GaussianRandomWalkProposalKernel::get_variables() const
{
  std::vector<std::string> variables;
  variables.reserve(this->proposal_info.size());
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    variables.push_back(i->first);
  }
  return variables;
}

GradientEstimatorOutput* GaussianRandomWalkProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> GaussianRandomWalkProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void GaussianRandomWalkProposalKernel::set_index(Index* index_in)
{
}

void GaussianRandomWalkProposalKernel::set_index_if_null(Index* index_in)
{
}

bool GaussianRandomWalkProposalKernel::can_be_evaluated() const
{
  return true;
}

void GaussianRandomWalkProposalKernel::set_data(Data* data_in)
{
  
}
