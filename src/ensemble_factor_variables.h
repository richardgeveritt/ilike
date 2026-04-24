#ifndef ENSEMBLEFACTORVARIABLES_H
#define ENSEMBLEFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{
  /**
   * @file ensemble_factor_variables.h
   * @brief Defines the Particle class.
   *
   * Provides particle functionality.
   *
   * @namespace ilike
   * @class Particle
   * @brief The particle class.
   */


class Particle;
class Index;
class EnsembleFactors;
class EnsembleKalmanInversion;

class EnsembleFactorVariables
{
  
public:
  
  /**
   * @brief Performs the ensemblefactorvariables operation.
   */
  EnsembleFactorVariables();
  //EnsembleFactorVariables(EnsembleFactors* ensemble_factors_in);
  /**
   * @brief Performs the ~ensemblefactorvariables operation.
   */
  virtual ~EnsembleFactorVariables();
  
  /**
   * @brief Performs the ensemblefactorvariables operation.
   *
   * @param another The Particle instance to copy from.
   */
  EnsembleFactorVariables(const EnsembleFactorVariables &another);
  
  /**
   * @brief Assignment operator for Particle.
   *
   * @param another The Particle instance to copy from.
   */
  void operator=(const EnsembleFactorVariables &another);
  /**
   * @brief Creates a deep copy of this Particle object.
   *
   * @return The result.
   */
  virtual EnsembleFactorVariables* duplicate() const=0;
  
  // return type is ok for all of the subclasses we have so far - one arma::mat per measurement
  /**
   * @brief Returns the measurement states for covariance.
   *
   * @return The result.
   */
  virtual std::vector<arma::rowvec> get_measurement_states_for_covariance() const=0;
  /**
   * @brief Returns the shifts.
   *
   * @param inverse_incremental_temperature The inverse incremental temperature.
   *
   * @return The result.
   */
  virtual std::vector<arma::colvec> get_shifts(double inverse_incremental_temperature) const=0;
  /**
   * @brief Returns the deterministic shifts.
   *
   * @return The result.
   */
  virtual std::vector<arma::colvec> get_deterministic_shifts() const=0;
  //virtual std::vector<arma::mat> get_unconditional_measurement_covariances(const std::vector<arma::mat> &Cyys,
  //                                                                         double inverse_incremental_temperature) const=0;
  
  virtual double evaluate_ensemble_likelihood_ratios(const Index* index,
                                                     double inverse_incremental_temperature,
                                                     const std::vector<arma::mat> &inv_sigma_precomps,
                                                     const std::vector<double> &log_det_precomps) const=0;
  
  virtual double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                               double inverse_incremental_temperature,
                                                               const std::vector<arma::mat> &inv_sigma_precomps,
                                                               const std::vector<double> &log_det_precomps) const=0;
  
  /*
   virtual double evaluate_ensemble_likelihood_ratios(const Index* index,
   double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters)=0;
   virtual double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
   double inverse_incremental_temperature,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual double evaluate_likelihoods(const Index* index) const=0;
  
  /**
   * @brief Performs the subsample evaluate likelihoods operation.
   *
   * @param index The index.
   *
   * @return The result.
   */
  virtual double subsample_evaluate_likelihoods(const Index* index) const=0;
  
  /*
   virtual double evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters)=0;
   virtual double subsample_evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters)=0;
   */
  
  virtual const EnsembleFactors* get_ensemble_factors() const=0;
  
  /**
   * @brief Sets the particle.
   *
   * @param particle_in The particle.
   */
  void set_particle(Particle* particle_in);
  /**
   * @brief Returns the particle.
   *
   * @return The result.
   */
  Particle* get_particle();
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  virtual void forget_you_were_already_written_to_file()=0;
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index) const=0;
  
  /**
   * @brief Closes any open file streams.
   */
  virtual void close_ofstreams()=0;
  
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
  friend EnsembleKalmanInversion;
  
  // not stored here
  /** @brief The particle. */
  Particle* particle;
  //EnsembleFactors* ensemble_factors;
  
  /**
   * @brief Copies the state of another Particle into this object.
   *
   * @param another The Particle instance to copy from.
   */
  void make_copy(const EnsembleFactorVariables &another);
  
};
}

#endif
