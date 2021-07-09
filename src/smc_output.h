#ifndef SMCOUTPUT_H
#define SMCOUTPUT_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <deque>

#include "likelihood_estimator_output.h"
#include "particles.h"

class SMC;

class SMCOutput : public LikelihoodEstimatorOutput
{

public:

  SMCOutput();
  virtual ~SMCOutput();

  SMCOutput(SMC* estimator_in,
            size_t lag_in,
            size_t lag_proposed_in);

  SMCOutput(const SMCOutput &another);
  void operator=(const SMCOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  SMCOutput* smc_duplicate() const;

  void continue_simulate(const Parameters &parameters);
  void estimate(const Parameters &parameters);

  void add_particles(const Particles &latest_particles);
  void add_proposed_particles(const Particles &latest_proposals);
  void add_weights(const arma::colvec &latest_unnormalised_log_weight_updates);

  void print(std::ostream &os) const;

protected:

  // Stored in ModelAndAlgorithm or in main.
  SMC* estimator;

  std::deque<Particles> all_particles;

  std::deque<Particles> all_proposed;

  std::deque<arma::colvec> unnormalised_log_weights;

  std::deque<arma::colvec> normalised_log_weights;

  std::vector<double> log_normalising_constant_ratios;

  size_t lag;

  size_t lag_proposed;

  void make_copy(const SMCOutput &another);

};

#endif
