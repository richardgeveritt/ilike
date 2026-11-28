#ifndef AVERAGEDLATENTSTEP_H
#define AVERAGEDLATENTSTEP_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <string>
#include <vector>
#include "ratio_step.h"
#include "ilike_header.h"

namespace ilike
{

class Index;

// MHAAR-RB as a state step of the ratio-estimating moves (dev/plan.md 3.22
// U5b, C312; Andrieu, Yildirim, Doucet & Chopin, arXiv 2101.01253, section
// 3.2 and Algorithm 4), for a model pi(theta, z) proportional to
// eta(theta) prod_t gamma_{t,theta}(z_t): every term reads at most one
// element z_t of the latent variable z. The default MHAAR-RB-0 (gamma_{t,theta,vartheta}
// = gamma_{t,theta}), with q_t the latent's own distribution (its statement)
// at the proposed parameters for the direction c = 1 and at the current ones
// for c = 2. With L the log target, Delta_{t,phi}(v) the change in L when
// z_t alone is set to v (at parameters phi), D_t(v) the same change in the
// factor that holds z's statement at q's parameters, and candidates
// v_t^(1) = z_t, v_t^(2..M) drawn from q, both directions give
//   log r = L(vartheta, z) - L(theta, z)
//           + sum_t [lse_i(Delta_{t,vartheta}(v_t^(i)) - D_t(v_t^(i)))
//                    - lse_i(Delta_{t,theta}(v_t^(i)) - D_t(v_t^(i)))]
// (the move adds the proposal's own terms), and the new latent z'_t =
// v_t^(k_t), k_t drawn with probabilities proportional to
// exp(Delta_{t,vartheta} - D_t). With refresh, a rejection sets z_t =
// v_t^(l_t), l_t with probabilities proportional to exp(Delta_{t,theta} - D_t)
// (the paper's optional stage, never rejected in the RB-0 case). The changes
// at each element are made in place (plates recompute only the terms reading
// it): O(MT) term evaluations at each of theta and vartheta.
class AveragedLatentStep : public RatioStep
{
public:

  // latent_factors: the factors reading z; own_factors: the factor(s)
  // holding z's statement (0-based); simulate: z drawn from its own
  // distribution given the other parameters
  AveragedLatentStep(const std::string &variable_in, size_t points_in, bool refresh_in,
                     SimulateIndependentProposalPtr simulate_in,
                     const std::vector<size_t> &latent_factors_in, const std::vector<size_t> &own_factors_in);
  AveragedLatentStep(const AveragedLatentStep &another);
  virtual ~AveragedLatentStep();
  RatioStep* duplicate() const;

  double draw(RandomNumberGenerator &rng, const Parameters &from, const Parameters &to) const;
  bool changes_state() const { return true; }
  bool uses_direction() const { return true; }
  bool refreshes() const { return this->refresh; }
  double draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed) const;
  double draw_state(RandomNumberGenerator &rng, const Index* index, const Particle &current, Particle &proposed,
                    int direction, Particle* refreshed, bool* refreshed_made) const;

protected:

  // the changes, at each element t and candidate i, of the factors reading
  // z (full) and of the factor holding its statement (own), at the
  // particle's parameters
  void element_changes(Particle &particle, const std::vector<arma::mat> &candidates, arma::mat &full, arma::mat &own) const;
  // z set to the chosen candidates, the factors reading it evaluated
  void install(Particle &particle, const std::vector<arma::mat> &candidates, const std::vector<size_t> &chosen) const;

  std::string variable;
  size_t points = 2;
  bool refresh = true;
  SimulateIndependentProposalPtr simulate = NULL;
  std::vector<size_t> latent_factors;
  std::vector<size_t> own_factors;
  std::shared_ptr<Index> latent_index;
  std::shared_ptr<Index> own_index;
};

}

#endif
