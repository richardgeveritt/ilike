#ifndef ANNEALEDPATH_H
#define ANNEALEDPATH_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>
#include <vector>
#include "distributions.h"
#include "parameters.h"

namespace ilike
{

// The functions of an unnormalised distribution gamma(x | theta)/Z(theta)
// over one variable (a field): its log density up to Z, one sweep of
// updates in place (and the same in reverse order), and an exact draw.
typedef double (*UnnormalisedPotentialPtr)(const arma::mat &field, const Parameters &parameters, const Data &data);
typedef void (*UnnormalisedSweepFunctionPtr)(RandomNumberGenerator &rng, Parameters &parameters, const Data &data);
typedef Parameters (*UnnormalisedExactFunctionPtr)(RandomNumberGenerator &rng, const Parameters &parameters, const Data &data);

// A path of distributions between two values of the parameters of a
// normalising constant, along which states are moved by MCMC kernels and
// weighted (dev/plan.md 3.22 U3c): the shared part of the annealed SAV/MAV
// estimate (C283), the estimate of log Z at a point (normalise, C289) and
// annealed exchange (C292). With the points
//   vartheta_k = start + k/(steps) (end - start), k = 0, ..., steps,
// (steps = bridges + 1) and states drawn from the distribution at start:
//  single(): one state, no resampling (annealed importance sampling): move
//   at each intermediate point, and return
//   sum_{k=0}^{steps-1} log gamma(u_k | vartheta_{k+1}) - log gamma(u_k | vartheta_k),
//   whose exponential is unbiased for Z(end)/Z(start);
//  population(): states moved together, with multinomial resampling after
//   each reweighting (SMC with MCMC moves), returning the sum over steps of
//   the log mean incremental weight, an unbiased estimate of log Z(end)/Z(start)
//   in the sense of the exponential.
// Engines: FieldPath (below: tight loops over one field's functions) and,
// later, the general SMC on the distribution's sub-model (U3c, with U4).
class AnnealedPath
{
public:
  virtual ~AnnealedPath() {}
  virtual AnnealedPath* duplicate() const=0;
  virtual double single(RandomNumberGenerator &rng, Parameters &state,
                        const Parameters &start, const Parameters &end) const=0;
  // threshold: resample when the effective sample size falls below
  // threshold x the number of states (1: at every step)
  virtual double population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                            const Parameters &start, const Parameters &end,
                            double threshold = 1.0) const=0;
};

// The fast engine: one field, its generated (or user-supplied) potential and
// sweeps. The kernel at a point is `sweeps` sweeps; with reversible, each is
// a sweep then a reverse sweep (a reversible kernel, which annealed exchange
// needs for exactness; importance sampling needs only invariance).
class FieldPath : public AnnealedPath
{
public:
  FieldPath();
  FieldPath(const std::string &variable_in,
            UnnormalisedPotentialPtr potential_in,
            UnnormalisedSweepFunctionPtr sweep_in,
            UnnormalisedSweepFunctionPtr reverse_sweep_in,
            const std::vector<std::string> &normaliser_in,
            size_t bridges_in,
            size_t sweeps_in,
            bool reversible_in,
            Data* data_in);
  AnnealedPath* duplicate() const;
  double single(RandomNumberGenerator &rng, Parameters &state,
                const Parameters &start, const Parameters &end) const;
  double population(RandomNumberGenerator &rng, std::vector<Parameters> &states,
                    const Parameters &start, const Parameters &end,
                    double threshold = 1.0) const;
  void set_data(Data* data_in) { this->data = data_in; }

protected:
  // the normaliser's parameters in `state` set to point k
  void set_point(Parameters &state, size_t k, const Parameters &start, const Parameters &end) const;
  void move(RandomNumberGenerator &rng, Parameters &state) const;

  std::string variable;
  UnnormalisedPotentialPtr potential = NULL;
  UnnormalisedSweepFunctionPtr sweep = NULL;
  UnnormalisedSweepFunctionPtr reverse_sweep = NULL;
  std::vector<std::string> normaliser;
  size_t bridges = 0;
  size_t sweeps = 1;
  bool reversible = false;
  // not stored here
  Data* data = NULL;
};

// log of the mean of exp(log_weights)
double log_mean_exp(const std::vector<double> &log_weights);

}

#endif
