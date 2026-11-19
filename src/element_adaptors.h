#ifndef ELEMENTADAPTORS_H
#define ELEMENTADAPTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "mcmc_adaptor.h"
#include "smc_adaptor.h"

namespace ilike
{

class ElementRandomWalkProposalKernel;

// Adaptation of moves on single elements (for (j in 1:J) move m =
// metropolis(norm_rw(theta[j], s), adapt = scale(…), smc_adapt = gaussian(c))):
// dev/plan.md, "adaptation of element moves". Each element's kernel has its
// own adaptor, so each element has its own step.

// Along the chain: the step adapted to the element's own acceptance rate,
// log s_{n+1} = log s_n + gamma_{n+1} (alpha_{n+1} - target), gamma_n =
// C (n+1)^-kappa (Andrieu & Thoms 2008, Algorithm 5, one element at a time);
// with averaging, the kernel uses the average of the log steps. The
// acceptance rate of a single-element move reflects the element's conditional
// spread, so this targets the conditional scale.
class ElementScaleAdaptor : public MCMCAdaptor
{
public:

  ElementScaleAdaptor();
  ElementScaleAdaptor(double target_in,
                      double gain_exponent_in,
                      double gain_constant_in,
                      bool averaging_in);
  virtual ~ElementScaleAdaptor();

  ElementScaleAdaptor(const ElementScaleAdaptor &another);
  void operator=(const ElementScaleAdaptor &another);
  MCMCAdaptor* duplicate() const;

  void adapt(const AdaptationStep &step);
  void append_state(std::vector<std::string> *names,
                    std::vector<double> *values,
                    bool final,
                    const std::string &prefix) const;

  bool watched(const std::string &name, std::vector<double> &values) const;
  std::vector<std::string> watchable() const;
  bool acceptance_target(double &target) const;
  std::string state_label() const;

protected:

  void make_copy(const ElementScaleAdaptor &another);
  ElementRandomWalkProposalKernel* kernel() const;

  double target;
  double gain_exponent;
  double gain_constant;
  bool averaging;
  size_t number_of_updates;
  double log_scale;
  double average_log_scale;
  bool started;
};

// From the SMC particles: the step is c times the element's conditional
// standard deviation given all the other parameters, 1/sqrt((Sigma^-1)_jj),
// with Sigma the particles' weighted covariance (shrunk towards its diagonal
// if singular; covariance_tools.h). c defaults to 2.38, the rule of thumb in
// one dimension. If it cannot be estimated, the step is left as it was.
class ElementSMCAdaptor : public SMCAdaptor
{
public:

  ElementSMCAdaptor();
  ElementSMCAdaptor(double constant_in);
  virtual ~ElementSMCAdaptor();

  ElementSMCAdaptor(const ElementSMCAdaptor &another);
  void operator=(const ElementSMCAdaptor &another);
  SMCAdaptor* duplicate() const;

  void smc_adapt(SMCOutput* current_state);
  void ensemble_adapt(EnsembleKalmanOutput* current_state);

protected:

  void make_copy(const ElementSMCAdaptor &another);
  // the step from the points (all the parameters; names and offsets as in
  // population_points) and their weights
  void adapt_from(const arma::mat &points,
                  const arma::colvec &normalised_weights,
                  const std::vector<std::string> &names,
                  const std::vector<size_t> &offsets,
                  size_t number_of_rows);

  double constant;
};

}

#endif
