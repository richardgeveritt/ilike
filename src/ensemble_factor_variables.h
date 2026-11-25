#ifndef ENSEMBLEFACTORVARIABLES_H
#define ENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{

class Particle;
class Index;
class EnsembleFactors;
class EnsembleKalmanInversion;

class EnsembleFactorVariables
{
  
public:
  
  EnsembleFactorVariables();
  // (C265) point these variables (and their measurement outputs) at
  // ensemble_factors_in, the copy owned by a copied filter; false if not possible
  virtual bool repoint_ensemble_factors(const EnsembleFactors* ensemble_factors_in) { return false; }
  // (C278) the member's predicted measurement for measurement factor j from
  // its current state (after the earlier factors' updates), without random
  // numbers; false if this kind of measurement cannot do it
  virtual bool recompute_measurement_state(size_t j, const Parameters &parameters) { return false; }
  virtual ~EnsembleFactorVariables();
  
  EnsembleFactorVariables(const EnsembleFactorVariables &another);
  
  void operator=(const EnsembleFactorVariables &another);
  virtual EnsembleFactorVariables* duplicate() const=0;
  
  // return type is ok for all of the subclasses we have so far - one arma::mat per measurement
  virtual std::vector<arma::rowvec> get_measurement_states_for_covariance() const=0;
  virtual std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const=0;
  virtual std::vector<arma::colvec> get_deterministic_shifts() const=0;
  
  virtual double evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature,
                                                     const std::vector<arma::mat> &inv_sigma_precomps,
                                                     const std::vector<double> &log_det_precomps) const=0;
  
  virtual double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature,
                                                               const std::vector<arma::mat> &inv_sigma_precomps,
                                                               const std::vector<double> &log_det_precomps) const=0;

  virtual double evaluate_likelihoods(const Index* index) const=0;
  
  virtual double subsample_evaluate_likelihoods(const Index* index) const=0;

  virtual const EnsembleFactors* get_ensemble_factors() const=0;
  
  void set_particle(Particle* particle_in);
  Particle* get_particle();
  
  virtual void forget_you_were_already_written_to_file()=0;
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index) const=0;
  
  virtual void close_ofstreams()=0;

protected:
  
  friend Particle;
  friend EnsembleKalmanInversion;
  
  // not stored here
  Particle* particle;
  
  void make_copy(const EnsembleFactorVariables &another);
  
};
}

#endif
