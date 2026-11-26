#ifndef INVOLUTIVEMCMC_H
#define INVOLUTIVEMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <atomic>
#include <string>
#include <vector>
#include "metropolis_hastings_mcmc.h"
#include "proposal_kernel.h"

namespace ilike
{

typedef Parameters (*InvolutionAuxiliaryPtr)(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data);
typedef double (*InvolutionLogDensityPtr)(const Parameters &parameters, const Data &data);
typedef Parameters (*InvolutionMapPtr)(const Parameters &parameters, const Data &data);

// The settings of a move written as an involution by the user (dev/plan.md
// 3.22 E2): the state variables it changes, the auxiliary variables' law
// (a draw and its log density), the map phi and the log of its Jacobian
// determinant, and how many proposals to check phi(phi(.)) = identity on.
// Not a proposal with a density: the recipe makes the move an InvolutiveMCMC,
// which builds its proposed particles through build() (ProposalKernel's
// proposed_from: the move's factors get new outputs, as for any proposal).
class InvolutionProposalKernel : public ProposalKernel
{
public:

  InvolutionProposalKernel();
  virtual ~InvolutionProposalKernel();
  InvolutionProposalKernel(const std::vector<std::string> &state_variables_in,
                           InvolutionAuxiliaryPtr auxiliary_in,
                           InvolutionLogDensityPtr log_auxiliary_in,
                           InvolutionMapPtr map_in,
                           InvolutionLogDensityPtr log_jacobian_in,
                           size_t check_in);
  InvolutionProposalKernel(const InvolutionProposalKernel &another);
  void operator=(const InvolutionProposalKernel &another);
  Kernel* duplicate() const;
  ProposalKernel* proposal_kernel_duplicate() const;

  // the proposed particle with the state variables in `state` (new outputs
  // for the factors in new_outputs_index)
  Particle build(const Particle &particle, const Parameters &state, const Index* new_outputs_index) const;

  void set_proposal_parameters(Parameters* proposal_parameters_in);
  GradientEstimatorOutput* simulate_gradient_estimator_output() const;
  std::vector<const ProposalKernel*> get_proposals() const;
  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  bool can_be_evaluated() const;
  void set_data(Data* data_in);

  std::vector<std::string> state_variables;
  InvolutionAuxiliaryPtr auxiliary = NULL;
  InvolutionLogDensityPtr log_auxiliary = NULL;
  InvolutionMapPtr map = NULL;
  InvolutionLogDensityPtr log_jacobian = NULL;   // (NULL: volume-preserving, or discrete)
  size_t check = 0;

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

  void make_copy(const InvolutionProposalKernel &another);
};

// A move written as an involution (Andrieu, Lee & Livingstone, arXiv
// 2012.14881; dev/plan.md 3.22 E2). With x the move's state variables, it
// draws auxiliary variables v ~ q(. | x), maps (x, v) to (x', v') = phi(x, v)
// with phi an involution (phi(phi(x, v)) = (x, v)), and accepts x' with
// probability min{1, r}, where
//   r = pi(x') q(v' | x') |det J_phi(x, v)| / [pi(x) q(v | x)],
// pi given by the move's factors (those that may read x). This is the
// Metropolis-Hastings update for the joint law pi(x) q(v | x) and the
// involution phi; v is not kept. The first `check` proposals check phi(phi(x,
// v)) = (x, v) (a wrong map is otherwise silent).
class InvolutiveMCMC : public MetropolisHastingsMCMC
{
public:

  InvolutiveMCMC();
  InvolutiveMCMC(size_t number_of_iterations_in, InvolutionProposalKernel* proposal_in, Data* data_in);
  InvolutiveMCMC(MCMCTermination* termination_in, InvolutionProposalKernel* proposal_in, Data* data_in);
  virtual ~InvolutiveMCMC();
  InvolutiveMCMC(const InvolutiveMCMC &another);
  void operator=(const InvolutiveMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  MetropolisHastingsMCMC* metropolis_hastings_mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng, const Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng, Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch) const;
  void adaptive_move_in_place(RandomNumberGenerator &rng, Particle &particle, Particle &scratch);
  Particle subsample_move(RandomNumberGenerator &rng, const Particle &particle) const;
  void diagnostic_names(std::vector<std::string> &names) const { names.clear(); }

protected:

  StandardMCMCOutput* initialise_mcmc_output() const;

  // the proposal and the decision: proposed_particle is set; returns whether accepted
  bool step(RandomNumberGenerator &rng, const Particle &particle, Particle &proposed_particle) const;

  void make_copy(const InvolutiveMCMC &another);

  // not stored here
  InvolutionProposalKernel* involution;
  Data* data;
  mutable std::atomic<size_t> checked;
};

}

#endif
