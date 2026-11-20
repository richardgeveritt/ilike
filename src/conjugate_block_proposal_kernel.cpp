#include "conjugate_block_proposal_kernel.h"

namespace ilike
{

// --- DirichletBlockProposalKernel -------------------------------------------

DirichletBlockProposalKernel::DirichletBlockProposalKernel()
:HMCProposalKernel()
{
}

DirichletBlockProposalKernel::~DirichletBlockProposalKernel()
{
}

DirichletBlockProposalKernel::DirichletBlockProposalKernel(const std::string &variable_in)
:HMCProposalKernel(variable_in, arma::mat(1, 1, arma::fill::ones), 1.0, 1)
{
}

DirichletBlockProposalKernel::DirichletBlockProposalKernel(const DirichletBlockProposalKernel &another)
:HMCProposalKernel(another)
{
}

void DirichletBlockProposalKernel::operator=(const DirichletBlockProposalKernel &another)
{
  if (this == &another)
    return;
  HMCProposalKernel::operator=(another);
}

Kernel* DirichletBlockProposalKernel::duplicate() const
{
  return new DirichletBlockProposalKernel(*this);
}

ProposalKernel* DirichletBlockProposalKernel::proposal_kernel_duplicate() const
{
  return new DirichletBlockProposalKernel(*this);
}

Parameters DirichletBlockProposalKernel::trajectory(RandomNumberGenerator &rng,
                                                    const Particle &particle,
                                                    double &log_ratio) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec p = current.get_colvec(this->variable);
  arma::uword K = p.n_elem;
  std::string name = stacked_display_name(this->variable);
  if (arma::any(p <= 0.0) || (std::abs(arma::accu(p) - 1.0) > 1e-8))
    Rcpp::stop("gibbs(" + name + "): the current value of " + name + " is not a probability vector with positive elements.");
  const double h = 0.5;
  double f0 = this->target_at(particle, current, p);
  arma::colvec D(K), L(K);
  for (arma::uword j = 0; j < K; ++j)
  {
    L[j] = std::log1p(p[j]*std::expm1(h));
    arma::colvec q = p;
    q[j] = q[j]*std::exp(h);
    q = q/arma::accu(q);
    D[j] = this->target_at(particle, current, q) - f0;
  }
  double S = arma::accu(D)/(h - arma::accu(L));
  arma::colvec c = (D + S*L)/h;
  arma::colvec alpha = c + 1.0;
  if (!alpha.is_finite() || arma::any(alpha <= 0.0))
    Rcpp::stop("gibbs(" + name + "): the log density of " + name + " given the rest is not that of a Dirichlet distribution.");
  arma::colvec x(K);
  for (arma::uword k = 0; k < K; ++k)
    x[k] = rgamma(rng, alpha[k], 1.0);
  x = x/arma::accu(x);
  // log q(p) - log q(p*), q = Dirichlet(alpha)
  log_ratio = arma::accu(c % (arma::log(p) - arma::log(x)));
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

// --- NormalInverseGammaBlockProposalKernel -----------------------------------

NormalInverseGammaBlockProposalKernel::NormalInverseGammaBlockProposalKernel()
:HMCProposalKernel(), mean_first(true)
{
}

NormalInverseGammaBlockProposalKernel::~NormalInverseGammaBlockProposalKernel()
{
}

NormalInverseGammaBlockProposalKernel::NormalInverseGammaBlockProposalKernel(const std::string &variable_in, bool mean_first_in)
:HMCProposalKernel(variable_in, arma::mat(1, 1, arma::fill::ones), 1.0, 1), mean_first(mean_first_in)
{
}

NormalInverseGammaBlockProposalKernel::NormalInverseGammaBlockProposalKernel(const NormalInverseGammaBlockProposalKernel &another)
:HMCProposalKernel(another), mean_first(another.mean_first)
{
}

void NormalInverseGammaBlockProposalKernel::operator=(const NormalInverseGammaBlockProposalKernel &another)
{
  if (this == &another)
    return;
  HMCProposalKernel::operator=(another);
  this->mean_first = another.mean_first;
}

Kernel* NormalInverseGammaBlockProposalKernel::duplicate() const
{
  return new NormalInverseGammaBlockProposalKernel(*this);
}

ProposalKernel* NormalInverseGammaBlockProposalKernel::proposal_kernel_duplicate() const
{
  return new NormalInverseGammaBlockProposalKernel(*this);
}

Parameters NormalInverseGammaBlockProposalKernel::trajectory(RandomNumberGenerator &rng,
                                                             const Particle &particle,
                                                             double &log_ratio) const
{
  const Parameters &current = particle.get_transformed_parameters(this);
  arma::colvec x0 = current.get_colvec(this->variable);
  if (x0.n_elem != 2)
    Rcpp::stop("gibbs([mean, variance]): the block has " + std::to_string(x0.n_elem) + " elements, not 2.");
  arma::uword im = this->mean_first ? 0 : 1, iv = this->mean_first ? 1 : 0;
  double mu0 = x0[im], v0 = x0[iv];
  std::string name = stacked_display_name(this->variable);
  if (!(v0 > 0.0))
    v0 = 1.0; // (a starting variance outside the support: fit around 1)
  // five points: (mu0, v0) and around it (the scale of mu: sqrt(v0))
  double s = std::sqrt(v0);
  const double mus[5] = {mu0, mu0 + s, mu0 - s, mu0, mu0 + s};
  const double vs[5] = {v0, v0, v0, v0*std::exp(0.5), v0*std::exp(-0.5)};
  arma::mat A(5, 5);
  arma::colvec f(5);
  for (int k = 0; k < 5; ++k)
  {
    arma::colvec x(2);
    x[im] = mus[k];
    x[iv] = vs[k];
    f[k] = this->target_at(particle, current, x);
    A(k, 0) = std::log(vs[k]);
    A(k, 1) = 1.0/vs[k];
    A(k, 2) = mus[k]/vs[k];
    A(k, 3) = mus[k]*mus[k]/vs[k];
    A(k, 4) = 1.0;
  }
  arma::colvec coefficients;
  bool solved = arma::solve(coefficients, A, f);
  double alpha = coefficients[0], beta = coefficients[1], gamma = coefficients[2], delta = coefficients[3];
  double shape = -alpha - 1.5, scale = -(beta - gamma*gamma/(4.0*delta));
  if (!solved || !(delta < 0.0) || !(shape > 0.0) || !(scale > 0.0) || !std::isfinite(shape) || !std::isfinite(scale))
    Rcpp::stop("gibbs(" + name + "): the log density given the rest is not that of a normal-inverse-gamma distribution.");
  double v = 1.0/rgamma(rng, shape, scale);
  double m = -gamma/(2.0*delta);
  double mu = rnorm(rng, m, std::sqrt(-v/(2.0*delta)));
  // log q(mu0, v0) - log q(mu, v), q proportional to the fitted density
  auto log_q = [&](double a, double b) {
    if (!(b > 0.0))
      return -arma::datum::inf;
    return alpha*std::log(b) + (beta + gamma*a + delta*a*a)/b;
  };
  log_ratio = log_q(x0[im], x0[iv]) - log_q(mu, v);
  arma::colvec x(2);
  x[im] = mu;
  x[iv] = v;
  Parameters output;
  output.set_stacked(this->variable, x, current);
  return output;
}

}
