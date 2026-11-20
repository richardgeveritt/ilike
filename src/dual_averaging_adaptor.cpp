#include "dual_averaging_adaptor.h"
#include "hmc_proposal_kernel.h"

namespace ilike
{

DualAveragingAdaptor::DualAveragingAdaptor()
:MCMCAdaptor(), target(0.8), gamma(0.05), t0(10.0), kappa(0.75), find_step(true)
{
}

DualAveragingAdaptor::DualAveragingAdaptor(double target_in, double gamma_in, double t0_in, double kappa_in, bool find_step_in,
                                           const std::string &mass_in)
:MCMCAdaptor(), target(target_in), gamma(gamma_in), t0(t0_in), kappa(kappa_in), find_step(find_step_in), mass(mass_in)
{
}

DualAveragingAdaptor::~DualAveragingAdaptor()
{
}

DualAveragingAdaptor::DualAveragingAdaptor(const DualAveragingAdaptor &another)
:MCMCAdaptor(another)
{
  this->make_copy(another);
}

void DualAveragingAdaptor::operator=(const DualAveragingAdaptor &another)
{
  if (this == &another)
    return;
  MCMCAdaptor::operator=(another);
  this->make_copy(another);
}

MCMCAdaptor* DualAveragingAdaptor::duplicate() const
{
  return new DualAveragingAdaptor(*this);
}

void DualAveragingAdaptor::make_copy(const DualAveragingAdaptor &another)
{
  this->target = another.target;
  this->gamma = another.gamma;
  this->t0 = another.t0;
  this->kappa = another.kappa;
  this->find_step = another.find_step;
  this->mass = another.mass;
  this->window_ends = another.window_ends;
  this->slow_start = another.slow_start;
  this->iteration = another.iteration;
  this->window_n = another.window_n;
  this->window_mean = another.window_mean;
  this->window_scatter = another.window_scatter;
  this->started = another.started;
  this->m = another.m;
  this->mu = another.mu;
  this->h_bar = another.h_bar;
  this->log_eps = another.log_eps;
  this->log_eps_bar = another.log_eps_bar;
}

HMCProposalKernel* DualAveragingAdaptor::kernel() const
{
  HMCProposalKernel* kernel = dynamic_cast<HMCProposalKernel*>(this->proposal);
  if (kernel == NULL)
    Rcpp::stop("dual_averaging adapts the step of ilike::hmc and ilike::nuts.");
  return kernel;
}

void DualAveragingAdaptor::set_step(double log_step)
{
  HMCProposalKernel* kernel = this->kernel();
  GaussianProposalInfo* info = kernel->get_gaussian_proposal_info(kernel->block());
  info->get_double_scale() = std::exp(2.0*log_step);
}

double DualAveragingAdaptor::find_reasonable_epsilon(RandomNumberGenerator &rng,
                                                     const HMCProposalKernel &kernel,
                                                     const Particle &particle,
                                                     double eps)
{
  const Parameters &base = particle.get_transformed_parameters(&kernel);
  arma::colvec x0 = base.get_colvec(kernel.block());
  arma::colvec r0 = kernel.draw_momentum(rng, x0.n_elem);
  arma::colvec g0 = kernel.gradient_at(particle, base, x0);
  double joint0 = kernel.target_at(particle, base, x0) - kernel.kinetic(r0);
  if (!std::isfinite(joint0) || !g0.is_finite())
    return eps;
  // log of the acceptance ratio of one leapfrog step of size e
  auto log_ratio = [&](double e) {
    arma::colvec r = r0 + (0.5*e)*g0;
    arma::colvec x = x0 + e*kernel.velocity(r);
    double L = kernel.target_at(particle, base, x);
    if (!std::isfinite(L))
      return -arma::datum::inf;
    arma::colvec g = kernel.gradient_at(particle, base, x);
    if (!g.is_finite())
      return -arma::datum::inf;
    r = r + (0.5*e)*g;
    double value = L - kernel.kinetic(r) - joint0;
    return std::isnan(value) ? -arma::datum::inf : value;
  };
  double lr = log_ratio(eps);
  // a = 1: double while the ratio is above 1/2; a = -1: halve while below
  double a = (lr > std::log(0.5)) ? 1.0 : -1.0;
  for (int k = 0; k < 100; ++k)
  {
    if (!(a*lr > -a*std::log(2.0)))
      break;
    eps = eps*std::pow(2.0, a);
    lr = log_ratio(eps);
  }
  return eps;
}

// Stan's schedule over a warm-up of W iterations (init_buffer 75, base
// window 25, term_buffer 50; scaled down to 15%, 75%, 10% when W < 150)
void DualAveragingAdaptor::plan_windows(size_t W)
{
  size_t init = 75, term = 50, base = 25;
  if (init + base + term > W)
  {
    init = size_t(0.15*double(W));
    term = size_t(0.1*double(W));
    base = W - init - term;
  }
  this->slow_start = init;
  this->window_ends.clear();
  size_t end_of_slow = W - term;
  size_t start = init, size = base;
  while (start < end_of_slow)
  {
    size_t end = start + size;
    // a next window that would not fit: this one runs to the end of the slow part
    if (end + 2*size > end_of_slow)
      end = end_of_slow;
    this->window_ends.push_back(end);
    start = end;
    size = 2*size;
  }
}

// a new start of dual averaging from the current step (or the step found
// from the particle's point)
void DualAveragingAdaptor::restart(RandomNumberGenerator* rng, const Particle* particle)
{
  HMCProposalKernel* kernel = this->kernel();
  double eps0 = kernel->step_size();
  if (this->find_step && (rng != NULL) && (particle != NULL))
    eps0 = find_reasonable_epsilon(*rng, *kernel, *particle, eps0);
  this->mu = std::log(10.0*eps0);
  this->log_eps = std::log(eps0);
  this->log_eps_bar = 0.0;
  this->h_bar = 0.0;
  this->m = 0;
  this->set_step(this->log_eps);
}

void DualAveragingAdaptor::adapt(const AdaptationStep &step)
{
  HMCProposalKernel* kernel = this->kernel();
  if (!this->started)
  {
    if (this->mass != "none")
    {
      size_t W = this->warm_up_length();
      if (W == 0)
        Rcpp::stop("dual_averaging(mass = " + this->mass + "): the warm-up's length must be known: until = iterations(n).");
      this->plan_windows(W);
    }
    this->restart(step.rng, step.previous);
    this->started = true;
  }
  this->iteration = this->iteration + 1;
  // the covariance in the slow windows (on the move's scale)
  if ((this->mass != "none") && (this->iteration > this->slow_start) && !this->window_ends.empty() &&
      (step.current != NULL))
  {
    arma::colvec z = step.current->get_transformed_parameters(kernel).get_colvec(kernel->block());
    if (this->window_n == 0)
    {
      this->window_mean = arma::zeros<arma::colvec>(z.n_elem);
      this->window_scatter = arma::zeros<arma::mat>(z.n_elem, z.n_elem);
    }
    this->window_n = this->window_n + 1;
    arma::colvec delta = z - this->window_mean;
    this->window_mean = this->window_mean + delta/double(this->window_n);
    this->window_scatter = this->window_scatter + delta*(z - this->window_mean).t();
    if (this->iteration == this->window_ends.front())
    {
      double n = double(this->window_n);
      if (this->window_n > 1)
      {
        arma::mat S = this->window_scatter/(n - 1.0);
        if (this->mass == "diagonal")
          S = arma::diagmat(S);
        arma::mat C = (n/(n + 5.0))*S + (1e-3*5.0/(n + 5.0))*arma::eye(S.n_rows, S.n_cols);
        kernel->get_gaussian_proposal_info(kernel->block())->set_covariance(C);
      }
      this->window_ends.erase(this->window_ends.begin());
      this->window_n = 0;
      this->restart(step.rng, step.current);
      return;
    }
  }
  this->m = this->m + 1;
  double md = double(this->m);
  double w = 1.0/(md + this->t0);
  this->h_bar = (1.0 - w)*this->h_bar + w*(this->target - step.acceptance_probability);
  this->log_eps = this->mu - std::sqrt(md)/this->gamma*this->h_bar;
  double eta = std::pow(md, -this->kappa);
  this->log_eps_bar = eta*this->log_eps + (1.0 - eta)*this->log_eps_bar;
  this->set_step(this->log_eps);
}

void DualAveragingAdaptor::on_stop()
{
  // the step from then on: the average of the iterates
  if (this->started)
  {
    this->log_eps = this->log_eps_bar;
    this->set_step(this->log_eps_bar);
  }
}

void DualAveragingAdaptor::append_state(std::vector<std::string> *names,
                                        std::vector<double> *values,
                                        bool final,
                                        const std::string &prefix) const
{
  std::string label = prefix + this->state_label();
  if (names != NULL)
  {
    names->push_back(label + ".step");
    names->push_back(label + ".step_average");
  }
  if (values != NULL)
  {
    values->push_back(this->started ? std::exp(this->log_eps) : this->kernel()->step_size());
    values->push_back(this->started ? std::exp(this->log_eps_bar) : this->kernel()->step_size());
  }
  // (C226) with mass windows, the final inverse mass matrix: its diagonal,
  // or its upper triangle (column by column) for dense
  if (final && (this->mass != "none"))
  {
    HMCProposalKernel* kernel = this->kernel();
    const arma::mat &C = kernel->get_gaussian_proposal_info(kernel->block())->get_covariance();
    for (arma::uword j = 0; j < C.n_cols; ++j)
      for (arma::uword i = 0; i <= j; ++i)
      {
        if ((this->mass == "diagonal") && (i != j))
          continue;
        if (names != NULL)
          names->push_back(label + ".mass[" + std::to_string(i + 1) + (this->mass == "dense" ? "," + std::to_string(j + 1) : std::string("")) + "]");
        if (values != NULL)
          values->push_back(C(i, j));
      }
  }
}

bool DualAveragingAdaptor::watched(const std::string &name, std::vector<double> &values) const
{
  // (the scale eps^2, as the scale adaptor of hmc and nuts; the average's)
  if (name != "scale")
    return false;
  values.assign(1, std::exp(2.0*this->log_eps_bar));
  return true;
}

std::vector<std::string> DualAveragingAdaptor::watchable() const
{
  return std::vector<std::string>(1, "scale");
}

bool DualAveragingAdaptor::acceptance_target(double &target_out) const
{
  target_out = this->target;
  return true;
}

std::string DualAveragingAdaptor::state_label() const
{
  return (this->proposal != NULL) ? stacked_display_name(this->kernel()->block()) : std::string("step");
}

}
