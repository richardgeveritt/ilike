#include "truncated_normal_independent_proposal_kernel.h"
#include "truncated_normal_distribution_factor.h"

namespace ilike
{

TruncatedNormalIndependentProposalKernel::TruncatedNormalIndependentProposalKernel()
:IndependentProposalKernel(), mean(0.0), sd(1.0), lower(-arma::datum::inf), upper(arma::datum::inf)
{
}

TruncatedNormalIndependentProposalKernel::~TruncatedNormalIndependentProposalKernel()
{
}

TruncatedNormalIndependentProposalKernel::TruncatedNormalIndependentProposalKernel(const std::string &variable_in, double mean_in,
                                                                                   double sd_in, double lower_in, double upper_in)
:IndependentProposalKernel(), variable(variable_in), mean(mean_in), sd(sd_in), lower(lower_in), upper(upper_in)
{
  if (!(sd_in > 0.0) || !(lower_in < upper_in))
    Rcpp::stop("ilike::tnorm: the standard deviation must be positive and lower < upper.");
}

TruncatedNormalIndependentProposalKernel::TruncatedNormalIndependentProposalKernel(const TruncatedNormalIndependentProposalKernel &another)
:IndependentProposalKernel(another)
{
  this->make_copy(another);
}

void TruncatedNormalIndependentProposalKernel::operator=(const TruncatedNormalIndependentProposalKernel &another)
{
  if (this == &another)
    return;
  IndependentProposalKernel::operator=(another);
  this->make_copy(another);
}

Kernel* TruncatedNormalIndependentProposalKernel::duplicate() const
{
  return new TruncatedNormalIndependentProposalKernel(*this);
}

ProposalKernel* TruncatedNormalIndependentProposalKernel::proposal_kernel_duplicate() const
{
  return new TruncatedNormalIndependentProposalKernel(*this);
}

IndependentProposalKernel* TruncatedNormalIndependentProposalKernel::independent_proposal_kernel_duplicate() const
{
  return new TruncatedNormalIndependentProposalKernel(*this);
}

void TruncatedNormalIndependentProposalKernel::make_copy(const TruncatedNormalIndependentProposalKernel &another)
{
  this->variable = another.variable;
  this->mean = another.mean;
  this->sd = another.sd;
  this->lower = another.lower;
  this->upper = another.upper;
}

arma::colvec TruncatedNormalIndependentProposalKernel::draw(RandomNumberGenerator &rng, arma::uword n) const
{
  arma::colvec x(n);
  for (arma::uword k = 0; k < n; ++k)
    x[k] = rtnorm(rng, this->mean, this->sd, this->lower, this->upper);
  return x;
}

double TruncatedNormalIndependentProposalKernel::evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  arma::colvec x = proposed_particle.get_colvec(this->variable);
  double output = 0.0;
  for (arma::uword k = 0; k < x.n_elem; ++k)
    output = output + log_truncated_normal(x[k], this->mean, this->sd, this->lower, this->upper);
  return output;
}

double TruncatedNormalIndependentProposalKernel::subsample_evaluate_independent_kernel(const Parameters &proposed_particle) const
{
  return this->evaluate_independent_kernel(proposed_particle);
}

Parameters TruncatedNormalIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng) const
{
  Parameters output;
  output[this->variable] = rtnorm(rng, this->mean, this->sd, this->lower, this->upper);
  return output;
}

Parameters TruncatedNormalIndependentProposalKernel::independent_simulate_like(RandomNumberGenerator &rng,
                                                                               const Parameters &shapes) const
{
  // one value, or each element of the block at the current point (C205)
  arma::uword n = 1;
  if (is_stacked_name(this->variable) || (shapes.find_vector(this->variable) != NULL))
    n = shapes.get_colvec(this->variable).n_elem;
  if (n == 1)
    return this->independent_simulate(rng);
  Parameters output;
  output.set_stacked(this->variable, this->draw(rng, n), shapes);
  return output;
}

Parameters TruncatedNormalIndependentProposalKernel::independent_simulate(RandomNumberGenerator &rng,
                                                                          const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng);
}

Parameters TruncatedNormalIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng) const
{
  return this->independent_simulate(rng);
}

Parameters TruncatedNormalIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                    const Parameters &conditioned_on_parameters) const
{
  return this->independent_simulate(rng);
}

Parameters TruncatedNormalIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                    const std::string &variable_in) const
{
  Parameters output;
  if (this->variable == variable_in)
    output[this->variable] = rtnorm(rng, this->mean, this->sd, this->lower, this->upper);
  return output;
}

Parameters TruncatedNormalIndependentProposalKernel::subsample_independent_simulate(RandomNumberGenerator &rng,
                                                                                    const std::string &variable_in,
                                                                                    const Parameters &conditioned_on_parameters) const
{
  return this->subsample_independent_simulate(rng, variable_in);
}

arma::mat TruncatedNormalIndependentProposalKernel::independent_gradient_of_log(const std::string &variable_in,
                                                                                const Parameters &proposed_particle)
{
  arma::mat x = proposed_particle[variable_in];
  if (variable_in != this->variable)
    return arma::zeros<arma::mat>(arma::size(x));
  return -(x - this->mean)/(this->sd*this->sd);
}

arma::mat TruncatedNormalIndependentProposalKernel::subsample_independent_gradient_of_log(const std::string &variable_in,
                                                                                          const Parameters &proposed_particle)
{
  return this->independent_gradient_of_log(variable_in, proposed_particle);
}

void TruncatedNormalIndependentProposalKernel::set_proposal_parameters(Parameters* proposal_parameters_in)
{
}

GradientEstimatorOutput* TruncatedNormalIndependentProposalKernel::simulate_gradient_estimator_output() const
{
  return NULL;
}

std::vector<const ProposalKernel*> TruncatedNormalIndependentProposalKernel::get_proposals() const
{
  std::vector<const ProposalKernel*> output;
  output.push_back(this);
  return output;
}

void TruncatedNormalIndependentProposalKernel::set_index(Index* index_in)
{
}

void TruncatedNormalIndependentProposalKernel::set_index_if_null(Index* index_in)
{
}

bool TruncatedNormalIndependentProposalKernel::can_be_evaluated() const
{
  return true;
}

void TruncatedNormalIndependentProposalKernel::set_data(Data* data_in)
{
}

}
