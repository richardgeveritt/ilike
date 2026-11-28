#ifndef INVOLUTIONSTEP_H
#define INVOLUTIONSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <atomic>
#include <memory>
#include <string>
#include <vector>
#include "ratio_step.h"
#include "ilike_header.h"
#include "annealed_path_kernel.h"

namespace ilike
{

typedef Parameters (*InvolutionAuxiliaryPtr)(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data);
typedef double (*InvolutionLogDensityPtr)(const Parameters &parameters, const Data &data);
typedef Parameters (*InvolutionMapPtr)(const Parameters &parameters, const Data &data);

// What an involution move draws and maps (dev/plan.md 3.22, C316): given xi
// (the state, with the proposed value of a model indicator as <k>_proposed
// when the move has one), the auxiliary variables u ~ q(. | x), log q(u | x)
// (with u in xi), the map (x, u) -> (x', u') (the state variables it
// changes, and u' under names that are not the model's), and log |det| of
// its Jacobian. The user's C++ functions (involution, rj), or built in
// (birth_death).
class InvolutionFunctions
{
public:
  virtual ~InvolutionFunctions() {}
  virtual InvolutionFunctions* duplicate() const=0;
  virtual Parameters draw(RandomNumberGenerator &rng, const Parameters &xi) const=0;
  virtual double log_density(const Parameters &xi) const=0;
  virtual Parameters map(const Parameters &xi) const=0;
  virtual double log_jacobian(const Parameters &xi) const { return 0.0; }
  // omega(k -> k', 1), the probability of the averaged form's direction c = 1
  virtual double first_direction_probability(double k, double k_new) const { return 0.5; }
  // whether the step checks that the map undoes itself (the user's maps)
  virtual bool checked() const { return true; }
  // whether the jump k -> k_new is one the functions can make
  virtual bool possible(double k, double k_new) const { return true; }
};

class UserInvolutionFunctions : public InvolutionFunctions
{
public:
  UserInvolutionFunctions(InvolutionAuxiliaryPtr auxiliary_in, InvolutionLogDensityPtr log_auxiliary_in,
                          InvolutionMapPtr map_in, InvolutionLogDensityPtr log_jacobian_in, const Data* data_in);
  InvolutionFunctions* duplicate() const { return new UserInvolutionFunctions(*this); }
  Parameters draw(RandomNumberGenerator &rng, const Parameters &xi) const;
  double log_density(const Parameters &xi) const;
  Parameters map(const Parameters &xi) const;
  double log_jacobian(const Parameters &xi) const;
protected:
  InvolutionAuxiliaryPtr auxiliary = NULL;
  InvolutionLogDensityPtr log_auxiliary = NULL;
  InvolutionMapPtr map_function = NULL;
  InvolutionLogDensityPtr log_jacobian_function = NULL; // (NULL: volume-preserving, or discrete)
  const Data* data = NULL;
};

// birth_death (dev/plan.md 3.22 U6): a birth inserts, at a position uniform
// among the k + 1, a new entry of each listed variable drawn from its
// statement given the other parameters (generated functions); a death removes
// a uniformly chosen entry; the later entries shift. u: the position and the
// new values (ilike_birth_position, ilike_birth_<v>), or the position removed
// (ilike_death_position). The positions' probabilities, 1/(k + 1) on both
// sides, cancel and are left out of the densities.
class BirthDeathFunctions : public InvolutionFunctions
{
public:
  BirthDeathFunctions(const std::string &indicator_in, size_t maximum_in, const std::vector<std::string> &variables_in,
                      SimulateIndependentProposalPtr simulate_in, EvaluateLogIndependentProposalPtr log_density_in);
  InvolutionFunctions* duplicate() const { return new BirthDeathFunctions(*this); }
  Parameters draw(RandomNumberGenerator &rng, const Parameters &xi) const;
  double log_density(const Parameters &xi) const;
  Parameters map(const Parameters &xi) const;
  double first_direction_probability(double k, double k_new) const { return (k_new > k) ? 1.0 : 0.0; }
  bool checked() const { return false; }
  bool possible(double k, double k_new) const;
protected:
  std::string indicator;
  size_t maximum = 1;
  std::vector<std::string> variables;
  SimulateIndependentProposalPtr simulate = NULL;
  EvaluateLogIndependentProposalPtr log_density_function = NULL;
};

// Every involution move as a state step of the ratio-estimating moves
// (dev/plan.md 3.22, C316; Andrieu, Lee & Livingstone, arXiv 2012.14881):
// the move's proposal supplies the part of the auxiliary variables chosen
// first (nothing for involution(…), whose proposal changes nothing; the model
// indicator's value k' for rj(…) and birth_death(…), dimension_jump), and the
// functions the rest: u ~ q(. | x), (x', u') = phi(x, u); the target part of
// the log acceptance ratio is
//   log pi(x') - log pi(x) + log q(u' | x') - log q(u | x) + log |J|.
// With an indicator, the move averages (points = N, RmJ-MCMC): its proposal
// fixes k' before the N draws.
class InvolutionStep : public RatioStep
{
public:

  // indicator: "" for involution(…); variables: the state variables the map changes
  // annealing_in: (stage B) with bridges, annealed reversible jump (Karagiannis
  // & Andrieu 2013): after the map, the state is moved along the path from
  // the start f0 to the end f1 (jump_factors.h) by the kernel's moves, and the
  // target part is the path's annealed importance sampling weight
  InvolutionStep(InvolutionFunctions* functions_in, const std::string &indicator_in, const std::vector<std::string> &variables_in, size_t check_in,
                 const AnnealedPathKernel &annealing_in = AnnealedPathKernel());
  InvolutionStep(const InvolutionStep &another);
  virtual ~InvolutionStep();
  RatioStep* duplicate() const;

  double draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const;
  bool changes_state() const { return true; }
  bool averages() const { return !this->indicator.empty(); }
  double first_direction_probability(const Particle &from, const Particle &to) const;
  double draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const;
  void link_kernel(const std::vector<const MCMC*> &cycle) { this->annealing.link_kernel(cycle, "rj"); }

protected:

  AnnealedPathKernel annealing;

  std::unique_ptr<InvolutionFunctions> functions;
  std::string indicator;
  std::vector<std::string> variables;
  size_t check = 10;
  std::shared_ptr<std::atomic<size_t>> checked_draws;
};

}

#endif
