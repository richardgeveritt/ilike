#ifndef GAUSSIANSMCADAPTOR_H
#define GAUSSIANSMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "smc_adaptor.h"
#include "gaussian_proposal_info.h"
#include "principal_components_info.h"

namespace ilike
{

// Adaptation of a Gaussian proposal from the SMC particles (dev/plan.md 3.2,
// B2): before each move, the proposal's covariance is set to the weighted
// sample covariance of the particles (for a principal-components proposal,
// its directions and variances to that covariance's eigenvectors and
// eigenvalues); for a mirror proposal its centre to their weighted mean; and,
// if a constant c is given, the scale to c^2/d (d the variable's dimension;
// c^2 for proposals that move one component or direction at a time). The
// weights are those of the particles for the current target, before
// resampling. If the covariance is not positive definite (e.g. too few
// distinct particles), the proposal keeps its previous covariance.
class GaussianSMCAdaptor : public SMCAdaptor
{

public:

  GaussianSMCAdaptor();

  // constant_in < 0: the scale is not changed
  GaussianSMCAdaptor(const std::string &variable_in,
                     double constant_in);

  virtual ~GaussianSMCAdaptor();

  GaussianSMCAdaptor(const GaussianSMCAdaptor &another);
  void operator=(const GaussianSMCAdaptor &another);
  SMCAdaptor* duplicate() const;

  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);

protected:

  void make_copy(const GaussianSMCAdaptor &another);

  // the weighted moments of the points (one row per point) into the kernel
  // points: all the parameters, the moved variable's d columns first
  void adapt_from(const arma::mat &points,
                  const arma::colvec &normalised_weights,
                  size_t d);

  std::string variable;
  double constant;

};

}

#endif
