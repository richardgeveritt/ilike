#include <iterator>
#include "uniform_random_walk_proposal_kernel.h"
#include "likelihood_estimator_output.h"
#include "likelihood_estimator.h"
#include "distributions.h"

namespace ilike
{

// the half-widths for a block of n elements (a variable or a stack): one per
// element, or a single one for every element (unif_rw on a vector or a list,
// C205)
static arma::colvec halfwidths_for(const std::string &variable, const arma::mat &halfwidth, arma::uword n)
{
  arma::colvec h = arma::vectorise(halfwidth);
  if (h.n_elem == n)
    return h;
  if (h.n_elem != 1)
    check_proposal_size(variable, n, h.n_elem, "half-width(s)", false);
  return arma::colvec(n, arma::fill::value(h[0]));
}
UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel()
:SymmetricProposalKernel()
{
}

UniformRandomWalkProposalKernel::~UniformRandomWalkProposalKernel()
{
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in)
:SymmetricProposalKernel()
{
  
  for (auto i=variable_names_in.begin();
       i!=variable_names_in.end();
       ++i)
  {
    this->proposal_info[*i] = arma::mat();
  }
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                 const double &halfwidth_in)
:SymmetricProposalKernel()
{
  
  this->proposal_info[variable_name_in] = halfwidth_in;
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                 const arma::colvec &halfwidth_in)
:SymmetricProposalKernel()
{
  
  this->proposal_info[variable_name_in] = halfwidth_in;
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::string &variable_name_in,
                                                                 const arma::mat &halfwidth_in)
:SymmetricProposalKernel()
{
  
  this->proposal_info[variable_name_in] = halfwidth_in;
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const std::vector<std::string> &variable_names_in,
                                                                 const std::vector<arma::mat> &halfwidths_in)
:SymmetricProposalKernel()
{
  
  for (size_t i=0;
       i<variable_names_in.size();
       ++i)
  {
    this->proposal_info[variable_names_in[i]] = halfwidths_in[i];
  }
}

UniformRandomWalkProposalKernel::UniformRandomWalkProposalKernel(const UniformRandomWalkProposalKernel &another)
:SymmetricProposalKernel(another)
{
  this->make_copy(another);
}

void UniformRandomWalkProposalKernel::operator=(const UniformRandomWalkProposalKernel &another)
{
  if(this == &another)
    return;
  
  SymmetricProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* UniformRandomWalkProposalKernel::duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

ProposalKernel* UniformRandomWalkProposalKernel::proposal_kernel_duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

SymmetricProposalKernel* UniformRandomWalkProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new UniformRandomWalkProposalKernel(*this));
}

void UniformRandomWalkProposalKernel::make_copy(const UniformRandomWalkProposalKernel &another)
{
  this->proposal_info = another.proposal_info;
}

double UniformRandomWalkProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                 const Particle &old_particle) const
{
  double output = 0.0;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    arma::colvec current_position = old_particle.get_transformed_parameters(this).get_colvec(i->first);
    arma::colvec h = halfwidths_for(i->first, i->second, current_position.n_elem);
    output = output + dunif(arma::mat(proposed_particle.get_transformed_parameters(this).get_colvec(i->first)),
                            arma::mat(current_position-h),
                            arma::mat(current_position+h));
  }
  return output;
}

double UniformRandomWalkProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                           const Particle &old_particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->specific_evaluate_kernel(proposed_particle, old_particle);
}

// simulate, written into the target (same random numbers and values)
void UniformRandomWalkProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                    const Particle &particle,
                                                    ProposalTarget &target) const
{
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    const Parameters &current = particle.get_transformed_parameters(this);
    arma::colvec current_position = current.get_colvec(i->first);
    arma::colvec h = halfwidths_for(i->first, i->second, current_position.n_elem);
    target.write(i->first, arma::vectorise(runif(rng, arma::mat(current_position-h), arma::mat(current_position+h))), current);
  }
}

Parameters UniformRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                     const Particle &particle) const
{
  Parameters output;
  for (auto i=this->proposal_info.begin();
       i!=this->proposal_info.end();
       ++i)
  {
    const Parameters &current = particle.get_transformed_parameters(this);
    arma::colvec current_position = current.get_colvec(i->first);
    arma::colvec h = halfwidths_for(i->first, i->second, current_position.n_elem);
    output.set_stacked(i->first, arma::vectorise(runif(rng, arma::mat(current_position-h), arma::mat(current_position+h))), current);
  }
  return output;
}

Parameters UniformRandomWalkProposalKernel::simulate(RandomNumberGenerator &rng,
                                                     const Particle &particle,
                                                     const Parameters &conditioned_on_parameters) const
{
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const Particle &particle,
                                                               const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->simulate(rng,particle);
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const std::string &variable,
                                                               const Particle &particle) const
{
  // no difference since size of data set does not impact on proposal
  auto found =  this->proposal_info.find(variable);
  
  Parameters output;
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec current_position = current.get_colvec(found->first);
  arma::colvec h = halfwidths_for(found->first, found->second, current_position.n_elem);
  output.set_stacked(found->first, arma::vectorise(runif(rng, arma::mat(current_position-h), arma::mat(current_position+h))), current);
  return output;
}

Parameters UniformRandomWalkProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                               const std::string &variable,
                                                               const Particle &particle,
                                                               const Parameters &conditioned_on_parameters) const
{
  // no difference since size of data set does not impact on proposal
  return this->subsample_simulate(rng,
                                  variable,
                                  particle);
}

arma::mat UniformRandomWalkProposalKernel::specific_gradient_of_log(const std::string &variable,
                                                                    const Particle &proposed_particle,
                                                                    const Particle &old_particle)
{
  Rcpp::stop("UniformRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

arma::mat UniformRandomWalkProposalKernel::specific_subsample_gradient_of_log(const std::string &variable,
                                                                              const Particle &proposed_particle,
                                                                              const Particle &old_particle)
{
  Rcpp::stop("UniformRandomWalkProposalKernel::specific_gradient_of_log - not written yet.");
}

void UniformRandomWalkProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
  
}

GradientEstimatorOutput* UniformRandomWalkProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> UniformRandomWalkProposalKernel::get_proposals() const
{
  std::vector<const  ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void UniformRandomWalkProposalKernel::set_index(Index* index_in)
{
}

void UniformRandomWalkProposalKernel::set_index_if_null(Index* index_in)
{
}

bool UniformRandomWalkProposalKernel::can_be_evaluated() const
{
  return true;
}

void UniformRandomWalkProposalKernel::set_data(Data* data_in)
{
  
}
}
