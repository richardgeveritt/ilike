#ifndef HMCPROPOSALKERNEL_H
#define HMCPROPOSALKERNEL_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <vector>

#include "proposal_kernel.h"
#include "particle.h"
#include "distributions.h"
#include "ilike_header.h"
#include "gaussian_proposal_info.h"

namespace ilike
{

class Index;

// Hamiltonian Monte Carlo (Neal 2011; Hoffman & Gelman, arXiv:1111.4246,
// Algorithm 1), as a proposal (dev/plan.md 3.18, C208). On the block x (a
// variable, or a stack of them), with C the "mass" argument (the inverse mass
// matrix M^-1, a covariance-scale matrix; 1 x 1 means c I) and step size
// eps (eps^2 is the GaussianProposalInfo's scale, so the existing adaptors
// adapt it): draw the momentum r ~ N(0, C^-1), make L leapfrog steps
//   r <- r + (eps/2) grad log pi(x);  x <- x + eps C r;  r <- r + (eps/2) grad log pi(x)
// with the gradients of the move's factors at each point, and propose x.
// The kinetic energy is K(r) = r' C r / 2. In a Metropolis-Hastings move the
// acceptance ratio is pi(x*)/pi(x) exp(K(r) - K(r*)): the proposed particle
// carries K(r) - K(r*) (its proposal_log_ratio) and a serial one more than
// the current particle's, and evaluate_kernel returns that difference as
// log q(old | new) and 0 as log q(new | old). In an unadjusted move the
// proposal is taken as it is (biased, as the unadjusted Langevin algorithm).
// With one step this is the Langevin (MALA) proposal with scale eps^2.
class HMCProposalKernel : public ProposalKernel
{

public:

  HMCProposalKernel();
  virtual ~HMCProposalKernel();

  HMCProposalKernel(const std::string &variable_in,
                    const arma::mat &covariance_in,
                    double step_in,
                    size_t steps_in);

  HMCProposalKernel(const HMCProposalKernel &another);
  void operator=(const HMCProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle,
                const Index* new_outputs_index) const;
  void move_into(RandomNumberGenerator &rng,
                 const Particle &particle,
                 const Index* new_outputs_index,
                 Particle &scratch) const;

  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

  GaussianProposalInfo* get_gaussian_proposal_info(const std::string &variable);

  // one trajectory from the particle: the proposed block's parameters, and
  // K(r) - K(r*) (minus infinity when the trajectory left the support)
  // (virtual: gibbs on a vector draws its normal conditional instead, C216)
  virtual Parameters trajectory(RandomNumberGenerator &rng,
                                const Particle &particle,
                                double &log_ratio) const;

  // The Hamiltonian system, for the No-U-Turn sampler (C210). On the
  // kernel's scale (transformed, with a transform): the block's name, the
  // step size, momenta r ~ N(0, C^-1), the velocity C r, K(r) = r' C r / 2.
  const std::string &block() const { return this->variable; }
  double step_size() const { return std::sqrt(this->info.get_double_scale()); }
  size_t get_max_depth() const { return this->max_depth; }
  size_t get_steps() const { return this->steps; }
  const Index* get_index() const { return this->index.get(); }
  void set_max_depth(size_t max_depth_in) { this->max_depth = max_depth_in; }
  // nuts's variant (C221): false, the paper's (slice sampling, the U-turn
  // criterion on the end points); true, Stan's (multinomial sampling, the
  // generalised criterion)
  bool get_stan_variant() const { return this->stan_variant; }
  void set_stan_variant(bool stan_variant_in) { this->stan_variant = stan_variant_in; }
  void check_size(arma::uword n) const;
  arma::colvec draw_momentum(RandomNumberGenerator &rng, arma::uword n) const;
  arma::colvec velocity(const arma::colvec &r) const;
  double kinetic(const arma::colvec &r) const;
  // the gradient of the move's log density at the block's value x (on the
  // kernel's scale: with a transform, x is z, and the gradient is that of
  // log pi(x(z)) + log|dx/dz|); `base` is the particle's parameters on that scale
  arma::colvec gradient_at(const Particle &particle, const Parameters &base, const arma::colvec &x) const;
  // the move's log density at x (its factors in `index`, plus log|dx/dz|
  // with a transform), with the particle at that point (its factors in
  // new_outputs_index given new outputs); minus infinity if x is not finite
  double log_density_at(const Particle &particle,
                        const Parameters &base,
                        const arma::colvec &x,
                        const Index* index,
                        const Index* new_outputs_index,
                        std::shared_ptr<Particle> &point) const;
  // the move's log density at x in the kernel's own index (log_density_at
  // without keeping the particle; for FindReasonableEpsilon, C219)
  double target_at(const Particle &particle, const Parameters &base, const arma::colvec &x) const;
  // log_density_at and gradient_at together, with each factor's value and
  // gradients in one pass where it can (C222): the particle at x in `point`,
  // the gradient on the kernel's scale in `gradient`
  double log_density_and_gradient_at(const Particle &particle,
                                     const Parameters &base,
                                     const arma::colvec &x,
                                     const Index* index,
                                     const Index* new_outputs_index,
                                     std::shared_ptr<Particle> &point,
                                     arma::colvec &gradient) const;
  // log|dx/dz| at the particle's point (0 without a transform)
  double log_jacobian(const Parameters &base) const;

protected:

  double specific_evaluate_kernel(const Particle &proposed_particle,
                                  const Particle &old_particle) const;
  double specific_subsample_evaluate_kernel(const Particle &proposed_particle,
                                            const Particle &old_particle) const;
  Parameters simulate(RandomNumberGenerator &rng,
                      const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const Particle &particle) const;
  Parameters subsample_simulate(RandomNumberGenerator &rng,
                                const std::string &variable,
                                const Particle &particle) const;
  arma::mat specific_gradient_of_log(const std::string &variable,
                                     const Particle &proposed_particle,
                                     const Particle &old_particle);
  arma::mat specific_subsample_gradient_of_log(const std::string &variable,
                                               const Particle &proposed_particle,
                                               const Particle &old_particle);

  void make_copy(const HMCProposalKernel &another);
  // the model-scale gradients of the block's variables, stacked and taken to
  // the kernel's scale at `point` (the parameters on that scale)
  arma::colvec stacked_gradient(const std::vector<arma::mat> &gradients, const Parameters &point, arma::uword n) const;

  std::string variable;
  GaussianProposalInfo info;
  size_t steps;
  size_t max_depth = 10; // (NUTS)
  bool stan_variant = false; // (NUTS, C221)
  std::unique_ptr<Index> index;
};

}

#endif
