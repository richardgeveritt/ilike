#include "normal_block_proposal_kernel.h"

namespace ilike
{

NormalBlockProposalKernel::NormalBlockProposalKernel()
:HMCProposalKernel()
{
}

NormalBlockProposalKernel::~NormalBlockProposalKernel()
{
}

NormalBlockProposalKernel::NormalBlockProposalKernel(const std::string &variable_in)
:HMCProposalKernel(variable_in, arma::mat(1, 1, arma::fill::ones), 1.0, 1)
{
}

NormalBlockProposalKernel::NormalBlockProposalKernel(const NormalBlockProposalKernel &another)
:HMCProposalKernel(another)
{
}

void NormalBlockProposalKernel::operator=(const NormalBlockProposalKernel &another)
{
  if (this == &another)
    return;
  HMCProposalKernel::operator=(another);
}

Kernel* NormalBlockProposalKernel::duplicate() const
{
  return new NormalBlockProposalKernel(*this);
}

ProposalKernel* NormalBlockProposalKernel::proposal_kernel_duplicate() const
{
  return new NormalBlockProposalKernel(*this);
}

bool NormalBlockProposalKernel::conditional(const Particle &particle,
                                            const Parameters &current,
                                            arma::colvec &mean,
                                            arma::mat &upper) const
{
  arma::colvec x0 = current.get_colvec(this->variable);
  arma::uword d = x0.n_elem;
  arma::colvec g0 = this->gradient_at(particle, current, x0);
  if (!g0.is_finite())
    return false;
  // h: 1% of the block's scale (any h is exact for a quadratic; this keeps
  // the differences well above rounding)
  double h = 1e-2*std::max(1.0, arma::abs(x0).max());
  arma::mat Q(d, d);
  for (arma::uword k = 0; k < d; ++k)
  {
    arma::colvec x = x0;
    x[k] = x[k] + h;
    arma::colvec g = this->gradient_at(particle, current, x);
    Q.col(k) = -(g - g0)/h;
  }
  if (!Q.is_finite())
    return false;
  Q = 0.5*(Q + Q.t());
  if (!arma::chol(upper, Q))
    return false;
  // mean = x0 + Q^-1 g0
  mean = x0 + arma::solve(arma::trimatu(upper), arma::solve(arma::trimatl(upper.t()), g0));
  return mean.is_finite();
}

Parameters NormalBlockProposalKernel::trajectory(RandomNumberGenerator &rng,
                                                 const Particle &particle,
                                                 double &log_ratio) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x0 = current.get_colvec(this->variable);
  arma::colvec mean;
  arma::mat upper;
  if (!this->conditional(particle, current, mean, upper))
    Rcpp::stop("gibbs(" + stacked_display_name(this->variable) + "): the log density of " + stacked_display_name(this->variable) +
               " given the rest is not a concave quadratic at its current value, so its full conditional is not normal "
               "(is every factor reading it normal in it?).");
  // x* = mean + U^-1 z, so that x* ~ N(mean, Q^-1)
  arma::colvec z = arma::vectorise(rnorm(rng, std::pair<size_t, size_t>(x0.n_elem, 1)));
  arma::colvec x = mean + arma::solve(arma::trimatu(upper), z);
  // log q(x0) - log q(x*), q = N(mean, Q^-1): -|U (x0 - mean)|^2/2 + |z|^2/2
  arma::colvec r0 = upper*(x0 - mean);
  log_ratio = -0.5*arma::dot(r0, r0) + 0.5*arma::dot(z, z);
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

}
