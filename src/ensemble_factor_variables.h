#ifndef ENSEMBLEFACTORVARIABLES_H
#define ENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

class Particle;
class Index;
class EnsembleFactors;
class IterativeEnsembleKalmanInversion;

class EnsembleFactorVariables
{

public:

  EnsembleFactorVariables();
  EnsembleFactorVariables(EnsembleFactors* ensemble_factors_in);
  virtual ~EnsembleFactorVariables();

  EnsembleFactorVariables(const EnsembleFactorVariables &another);

  void operator=(const EnsembleFactorVariables &another);
  virtual EnsembleFactorVariables* duplicate() const=0;
  
  // return type is ok for all of the subclasses we have so far - one arma::mat per measurement
  virtual std::vector<arma::rowvec> get_measurement_states_for_covariance() const=0;
  virtual std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const=0;
  virtual std::vector<arma::colvec> get_deterministic_shifts() const=0;
  virtual std::vector<arma::mat> get_kalman_gains(const std::vector<arma::mat> &Cxys,
                                                  const std::vector<arma::mat> &Cyys,
                                                  double inverse_incremental_temperature) const=0;
  virtual std::vector<arma::colvec*> get_measurements() const=0;
  
  virtual std::vector<arma::mat> get_adjustments(const arma::mat &Zf,
                                                 const arma::mat &Ginv,
                                                 const arma::mat &Ftranspose,
                                                 const std::vector<arma::mat> &Vs,
                                                 double inverse_incremental_temperature) const=0;
  
  virtual double evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature)=0;
  virtual double evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature,
                                             const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature,
                                                       const Parameters &conditioned_on_parameters)=0;
  
  virtual double evaluate_likelihoods(const Index* index)=0;
  virtual double evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)=0;
  
  EnsembleFactors* get_ensemble_factors();
  
  void set_particle(Particle* particle_in);
  Particle* get_particle();
  
  //void set_temperature(double temperature_in);
  
  /*
  virtual void evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  virtual void evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                     const Parameters &conditioned_on_parameters)=0;
  virtual double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  virtual double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                      const Parameters &conditioned_on_parameters)=0;
  virtual double evaluate_likelihoods(const Index* index)=0;
  virtual double evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)=0;
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_likelihoods(const Index* index)=0;
  virtual double subsample_evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat direct_get_gradient_of_log(const std::string &variable)=0;
  virtual arma::mat direct_get_gradient_of_log(const std::string &variable,
                                               const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const std::string &variable)=0;
  virtual arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable)=0;
  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable,
                                               const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable)=0;
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)=0;
  */

protected:
  
  friend Particle;
  friend IterativeEnsembleKalmanInversion;
  
  // not stored here
  Particle* particle;
  EnsembleFactors* ensemble_factors;

  void make_copy(const EnsembleFactorVariables &another);

};

#endif
