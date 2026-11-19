#ifndef GAUSSIANMCMCADAPTOR_H
#define GAUSSIANMCMCADAPTOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <string>

#include "mcmc_adaptor.h"
#include "gaussian_proposal_info.h"
#include "principal_components_info.h"

namespace ilike
{

// Adaptation of a Gaussian proposal, following Andrieu & Thoms (2008),
// Section 5 (dev/plan.md 3.2). With gain gamma_n = C (n+1)^(-kappa) at the
// n-th update (n = 1, 2, ...), and mu, Sigma the running estimates of the
// target's mean and covariance:
//   AM        (Algorithm 2): mu, Sigma from the chain; scale fixed at 2.38^2/d
//             (for a mirror proposal: the scale given).
//   RB_AM     (Algorithm 3): as AM, Rao-Blackwellised: X_{i+1} replaced by its
//             expectation over the accept/reject step,
//             alpha Y_{i+1} + (1 - alpha) X_i.
//   SCALED_AM (Algorithm 4): mu, Sigma as AM and the scale by Robbins-Monro,
//             log lambda += gamma (alpha - target).
//   SCALE     the scale recursion of Algorithm 4 alone (Sigma fixed).
//   COMPONENTWISE_AM (Algorithm 5, with a componentwise kernel): the moved
//             component k is the one where Y_{i+1} differs from X_i;
//             log lambda_k += gamma (alpha - target); mu, Sigma as AM.
//   PC_EIGEN, PC_ONLINE (Algorithm 8, with a principal-components kernel):
//             the direction l moved is the one along which Y_{i+1} - X_i
//             lies; log ell_l += gamma (alpha - target); mu, Sigma as AM; the
//             directions and the variances along them from Sigma, by an
//             eigen-decomposition (PC_EIGEN) or by one step of the online PCA
//             recursion of Section 5.3.2 (PC_ONLINE).
// The proposal's covariance is lambda Sigma: the scale goes into the kernel
// every iteration; Sigma (which needs a Cholesky factorisation) or the
// principal directions every refresh_every iterations (0: every d iterations,
// d the variable's dimension; every iteration for PC_ONLINE). For a mirror
// proposal the mean mu goes into the kernel every iteration too.
// With averaging (Polyak & Juditsky 1992), the recursions above run as they
// are, but the kernel uses the averages of their iterates, theta-bar_n =
// (theta_0 + ... + theta_n)/(n+1) (scales averaged on the log scale, the
// scale of their recursion): with a gain decreasing more slowly than 1/n
// (kappa < 1), the average converges at the optimal rate.
enum class GaussianAdaptationMethod { AM, RB_AM, SCALED_AM, SCALE, COMPONENTWISE_AM, PC_EIGEN, PC_ONLINE };

// The mirror proposal's default target acceptance rate for a variable of
// dimension d: the rate that maximised the effective sample size of the
// worst functional (x_j^2; ESS of x_j is much larger, the move being
// antithetic) on Gaussian targets with the exact centre and covariance, in
// experiments for d = 1, ..., 100 (dev/plan.md 3.2, C171); linear in log d
// between the dimensions tried, 0.23 beyond.
double mirror_default_target_acceptance(size_t d);

class GaussianMCMCAdaptor : public MCMCAdaptor
{

public:

  GaussianMCMCAdaptor();

  // target_acceptance_in 0: the mirror proposal's default for the
  // variable's dimension (mirror_default_target_acceptance)
  GaussianMCMCAdaptor(const std::vector<std::string> &variable_names_in,
                      GaussianAdaptationMethod method_in,
                      double target_acceptance_in,
                      double gain_exponent_in,
                      double gain_constant_in,
                      size_t refresh_every_in,
                      bool averaging_in = false);

  virtual ~GaussianMCMCAdaptor();

  GaussianMCMCAdaptor(const GaussianMCMCAdaptor &another);
  void operator=(const GaussianMCMCAdaptor &another);
  MCMCAdaptor* duplicate() const;

  void adapt(const AdaptationStep &step);

  void append_state(std::vector<std::string> *names,
                    std::vector<double> *values,
                    bool final,
                    const std::string &prefix) const;

  bool watched(const std::string &name, std::vector<double> &values) const;
  std::vector<std::string> watchable() const;
  std::vector<std::string> default_watch() const;
  bool acceptance_target(double &target) const;
  std::string state_label() const;

protected:

  void make_copy(const GaussianMCMCAdaptor &another);

  // looks up the kernel's parameters for each variable
  void bind();

  struct VariableState
  {
    std::string variable;
    GaussianProposalInfo* info; // not owned: inside the kernel adapted
    arma::colvec* component_scales; // not owned: componentwise kernels only
    PrincipalComponentsInfo* components; // not owned: principal-components kernels only
    arma::colvec mean;
    arma::mat covariance;
    double log_scale;
    arma::colvec log_component_scales; // componentwise scales, or the directions' scales
    arma::mat pca_loadings; // W of the online PCA recursion
    double target; // the target acceptance rate for this variable
    // the averages of the iterates (with averaging)
    double average_log_scale;
    arma::colvec average_mean;
    arma::mat average_covariance;
    arma::colvec average_log_component_scales;
    bool started;
    // conditional scaling of a move on part of the parameters (dev/plan.md):
    // the other non-fixed vector parameters, and the running statistics of
    // z = [x; r] (x this variable on the kernel's scale, r the rest): with
    // implicit weights w_k of the recursion, A = sum w_k y_k, C = sum w_k^2
    // y_k, B = sum w_k^2 y_k o y_k (y_k the outer products of the deviations),
    // W = sum w_k, Q = sum w_k^2 (for the weighted Ledoit-Wolf intensity)
    std::vector<std::string> rest;
    arma::colvec joint_mean;
    arma::mat joint_sum;
    arma::mat joint_weighted_products;
    arma::mat joint_squared_products;
    double weight_sum = 0.0;
    double squared_weight_sum = 0.0;
    arma::mat average_joint_covariance;
    arma::mat conditional_covariance; // the last one formed (for the output)
  };

  // the principal directions of the current Sigma, into the kernel
  void refresh_principal_components(VariableState &state,
                                    const arma::mat &covariance);

  bool uses_global_scale() const;
  // z = [x; r] at a point
  arma::colvec joint_point(const VariableState &state, const Particle &particle) const;
  // the conditional covariance of x given the rest (shrunk if singular), the
  // regression matrix Sigma_xr Sigma_rr^-1 over all of r (zero for the
  // coordinates of r that do not vary), and the means; FALSE if not yet estimable
  bool conditional_estimate(const VariableState &state,
                            arma::mat &conditional,
                            arma::mat &regression,
                            arma::colvec &x_mean,
                            arma::colvec &rest_mean) const;
  bool estimates_moments() const;

  std::vector<VariableState> states;
  GaussianAdaptationMethod method;
  double target_acceptance;
  double gain_exponent;
  double gain_constant;
  size_t refresh_every;
  size_t number_of_updates;
  bool writes_mean; // mirror proposals
  bool averaging; // Polyak-Juditsky averaging
  // the kernel the pointers were looked up in
  const ProposalKernel* bound_to;

};

}

#endif
