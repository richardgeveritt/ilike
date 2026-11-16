#include <cmath>
#include "principal_components_proposal_kernel.h"
#include "distributions.h"

namespace ilike
{

PrincipalComponentsProposalKernel::PrincipalComponentsProposalKernel()
:GaussianRandomWalkProposalKernel()
{
}

PrincipalComponentsProposalKernel::~PrincipalComponentsProposalKernel()
{
}

PrincipalComponentsProposalKernel::PrincipalComponentsProposalKernel(const std::string &variable_name_in,
                                                                     const arma::mat &covariance_in,
                                                                     double scale_in,
                                                                     size_t number_of_directions)
:GaussianRandomWalkProposalKernel(variable_name_in, covariance_in, scale_in)
{
  this->variable = variable_name_in;
  size_t d = covariance_in.n_rows;
  size_t m = number_of_directions;
  if (m==0)
    m = d;
  if (m>d)
    Rcpp::stop("ilike::pc_rw: more directions than the dimension of " + variable_name_in + ".");
  arma::colvec values;
  arma::mat vectors;
  arma::eig_sym(values, vectors, (covariance_in + covariance_in.t())/2.0);
  // eig_sym sorts ascending: the first m principal directions are the last m
  this->components.directions = arma::fliplr(vectors.tail_cols(m));
  this->components.eigenvalues = arma::flipud(values.tail(m));
  this->components.scales = arma::colvec(m, arma::fill::value(scale_in));
}

PrincipalComponentsProposalKernel::PrincipalComponentsProposalKernel(const PrincipalComponentsProposalKernel &another)
:GaussianRandomWalkProposalKernel(another)
{
  this->make_copy(another);
}

void PrincipalComponentsProposalKernel::operator=(const PrincipalComponentsProposalKernel &another)
{
  if(this == &another)
    return;
  GaussianRandomWalkProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* PrincipalComponentsProposalKernel::duplicate() const
{
  return( new PrincipalComponentsProposalKernel(*this));
}

ProposalKernel* PrincipalComponentsProposalKernel::proposal_kernel_duplicate() const
{
  return( new PrincipalComponentsProposalKernel(*this));
}

SymmetricProposalKernel* PrincipalComponentsProposalKernel::symmetric_proposal_kernel_duplicate() const
{
  return( new PrincipalComponentsProposalKernel(*this));
}

void PrincipalComponentsProposalKernel::make_copy(const PrincipalComponentsProposalKernel &another)
{
  this->variable = another.variable;
  this->components = another.components;
}

PrincipalComponentsInfo* PrincipalComponentsProposalKernel::get_principal_components(const std::string &variable_in)
{
  if (variable_in!=this->variable)
    return NULL;
  return &this->components;
}

arma::colvec PrincipalComponentsProposalKernel::move_along_one_direction(RandomNumberGenerator &rng,
                                                                         const Particle &particle) const
{
  arma::colvec point = particle.get_transformed_parameters(this).get_colvec(this->variable);
  if (point.n_elem!=this->components.directions.n_rows)
    Rcpp::stop("ilike::pc_rw: the covariance does not match the dimension of " + this->variable + ".");
  size_t m = this->components.directions.n_cols;
  // l uniform on 0, ..., m-1, then one normal increment along w_l
  size_t l = size_t(std::floor(runif(rng)*double(m)));
  if (l>=m)
    l = m-1;
  double sd = std::sqrt(this->components.scales[l]*this->components.eigenvalues[l]);
  point = point + (sd*rnorm(rng))*this->components.directions.col(l);
  return point;
}

Parameters PrincipalComponentsProposalKernel::simulate(RandomNumberGenerator &rng,
                                                       const Particle &particle) const
{
  Parameters output;
  output[this->variable] = this->move_along_one_direction(rng, particle);
  return output;
}

void PrincipalComponentsProposalKernel::simulate_into(RandomNumberGenerator &rng,
                                                      const Particle &particle,
                                                      ProposalTarget &target) const
{
  target.matrix(this->variable) = this->move_along_one_direction(rng, particle);
}

Parameters PrincipalComponentsProposalKernel::subsample_simulate(RandomNumberGenerator &rng,
                                                                 const Particle &particle) const
{
  return this->simulate(rng, particle);
}

double PrincipalComponentsProposalKernel::specific_evaluate_kernel(const Particle &proposed_particle,
                                                                   const Particle &old_particle) const
{
  Rcpp::stop("ilike::pc_rw: its density is singular (one direction moves); use it in an m_proposal, where it cancels.");
}

double PrincipalComponentsProposalKernel::specific_evaluate_kernel_from_parameters(const Parameters &proposed_parameters,
                                                                                   const Parameters &old_parameters) const
{
  Rcpp::stop("ilike::pc_rw: its density is singular (one direction moves); use it in an m_proposal, where it cancels.");
}

double PrincipalComponentsProposalKernel::specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                                                             const Particle &old_particle) const
{
  Rcpp::stop("ilike::pc_rw: its density is singular (one direction moves); use it in an m_proposal, where it cancels.");
}

}
