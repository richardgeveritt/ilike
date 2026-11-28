#ifndef INPLACEELEMENTSCANMCMC_H
#define INPLACEELEMENTSCANMCMC_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <string>
#include <vector>
#include "deterministic_scan_mcmc.h"
#include "element_random_walk_proposal_kernel.h"

namespace ilike
{

class Index;

// A scan of Metropolis random-walk moves on the elements of one variable,
// made in place (dev/plan.md I4c-2). For each element: the element is changed
// in the current particle; each plate factor in the move's index recomputes
// only the terms reading that element (keeping their old values), and every
// other factor in the index gets a new output (keeping the old one); then
// the proposal is accepted, or the element, the terms and the outputs are
// restored. No particle is copied, so an update costs the terms it touches
// (and O(log n) to re-add each plate's terms). The random numbers, and so
// the chain, are those of the same moves made through proposed particles:
// those moves (the element moves this scan holds, as a DeterministicScanMCMC)
// are used instead wherever the in-place path does not apply (factors not
// held as a vector of outputs, ensemble factors, subsampling).
// With `states` given (gibbs(s[j]), I4c-3) each element is instead drawn
// from its full conditional over those states: each state is set in turn
// and the factors in the move's index evaluated (each costing the terms
// reading the element), then a state is drawn and kept. The factors must be
// exact (a Gibbs draw needs the exact conditional); there is no copy-based
// equivalent, so the in-place path must apply.
// With step normal_conditional (gibbs on a normal mean, C206) the element's
// full conditional is normal, so its log density is exactly quadratic: it is
// evaluated (as for states) at three values, the quadratic gives the
// conditional's mean and variance, and a draw from it is accepted with the
// Metropolis-Hastings ratio of that independent proposal (1 up to rounding;
// computed, so the chain stays exact whatever the rounding).
// With step gamma_, invgamma_ or beta_conditional (C215) the log full
// conditional is exactly a phi1(v) + b phi2(v) + c, with (phi1, phi2) =
// (log v, v), (log v, 1/v) or (log v, log(1 - v)): evaluated at three values
// inside the support, it gives a and b, hence the gamma, inverse gamma or
// beta conditional, drawn and accepted in the same way.
class InPlaceElementScanMCMC : public DeterministicScanMCMC
{
public:

  InPlaceElementScanMCMC();
  // element_moves_in: one MetropolisMCMC with an ElementRandomWalkProposalKernel
  // per element, in the same order (takes ownership)
  InPlaceElementScanMCMC(MCMCTermination* termination_in,
                         const std::vector<MCMC*> &element_moves_in,
                         const std::string &variable_in,
                         const std::vector<size_t> &rows_in,
                         const std::vector<size_t> &columns_in,
                         const std::vector<double> &scales_in,
                         ElementStep step_in,
                         const std::vector<double> &states_in = std::vector<double>());
  virtual ~InPlaceElementScanMCMC();

  InPlaceElementScanMCMC(const InPlaceElementScanMCMC &another);
  void operator=(const InPlaceElementScanMCMC &another);
  Kernel* duplicate() const;
  MCMC* mcmc_duplicate() const;
  DeterministicScanMCMC* deterministic_scan_mcmc_duplicate() const;

  Particle move(RandomNumberGenerator &rng,
                const Particle &particle) const;
  // the scan reads the factors' cached values in every mode (a Metropolis
  // denominator, a Gibbs conditional); a gibbs scan holds no element moves,
  // so the scan's default (asking its moves) would say it does not, and an
  // SMC would then not reset the particle's target from the cache after the
  // move (C206)
  bool reads_cached_factors() const { return true; }
  void move_in_place(RandomNumberGenerator &rng,
                     Particle &particle) const;
  void move_in_place(RandomNumberGenerator &rng,
                     Particle &particle,
                     Particle &scratch) const;

  void set_index(Index* index_in);
  void set_index_if_null(Index* index_in);
  void set_cut_index(const Index* cut_index_in);
  // (U5a) the elements in reverse order; each element's update (a Metropolis
  // move or a draw from its full conditional) is reversible
  void move_in_place_reversed(RandomNumberGenerator &rng, Particle &particle) const;
  bool reversible_parts(std::string &why) const { return true; }

  // adaptation (dev/plan.md, "adaptation of element moves"): along the chain
  // each element's kernel adapts its step from the element's acceptance
  // probability (an ElementScaleAdaptor); from SMC particles every element's
  // step is set at once to c times its conditional sd given all the other
  // parameters (one precision matrix for the scan)
  void adaptive_move_in_place(RandomNumberGenerator &rng,
                              Particle &particle,
                              Particle &scratch);
  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);
  void set_smc_constant(double constant_in) { this->smc_constant = constant_in; }
  void append_adaptation_state(std::vector<std::string> *names,
                               std::vector<double> *values,
                               bool final,
                               const std::string &prefix) const;

protected:

  // false if the in-place path does not apply to this particle; with adapt,
  // each element's kernel adapts after its update
  bool scan_in_place(RandomNumberGenerator &rng,
                     Particle &particle,
                     bool adapt = false) const;
  // element e's kernel (NULL for gibbs scans), and its current step
  ElementRandomWalkProposalKernel* element_kernel(size_t e) const;
  double step_scale(size_t e) const;
  // every element's step from the points of a population (all parameters)
  void set_steps_from(const arma::mat &points,
                      const arma::colvec &normalised_weights,
                      const std::vector<std::string> &names,
                      const std::vector<size_t> &offsets,
                      size_t number_of_rows);
  double smc_constant = 0.0; // > 0: smc_adapt = gaussian(c) on the elements
  mutable bool reversed_order = false; // (U5a: move_in_place_reversed)

  std::string variable;
  std::vector<size_t> rows;
  std::vector<size_t> columns;
  std::vector<double> scales;
  ElementStep step_kind;
  std::vector<double> states;
  std::unique_ptr<Index> index;
  // a Gibbs draw (over states, or from a normal full conditional)
  bool is_gibbs() const
  {
    return !this->states.empty() || (this->step_kind == ElementStep::normal_conditional) ||
           (this->step_kind == ElementStep::gamma_conditional) || (this->step_kind == ElementStep::invgamma_conditional) ||
           (this->step_kind == ElementStep::beta_conditional);
  }

  void make_copy(const InPlaceElementScanMCMC &another);
};

}

#endif
