#ifndef FACTORVARIABLES_H
#define FACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

namespace ilike
{
  /**
   * @file factor_variables.h
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
class Factors;

class FactorVariables
{

public:

  /**
   * @brief Performs the factorvariables operation.
   */
  FactorVariables();
  /**
   * @brief Performs the ~factorvariables operation.
   */
  virtual ~FactorVariables();
  
  /**
   * @brief Performs the factorvariables operation.
   *
   * @param particle_in The particle.
   */
  FactorVariables(Particle* particle_in);

  /**
   * @brief Performs the factorvariables operation.
   *
   * @param another The Particle instance to copy from.
   */
  FactorVariables(const FactorVariables &another);

  /**
   * @brief Assignment operator for Particle.
   *
   * @param another The Particle instance to copy from.
   */
  void operator=(const FactorVariables &another);
  /**
   * @brief Creates a deep copy of this Particle object.
   *
   * @return The result.
   */
  virtual FactorVariables* duplicate() const=0;

  // Factor variables for a proposed point new_parameters: the factors in
  // new_outputs_index get new outputs (initialised and simulated at
  // new_parameters, not yet evaluated); every other factor's output is a copy
  // of this object's, keeping its cached value (and, for a random factor, its
  // estimate). Used by MCMC moves, whose other factors either cancel in the
  // acceptance ratio or are refreshed after acceptance. The caller owns the
  // result and must call set_particle on it.
  virtual FactorVariables* copy_with_new_outputs(const Parameters &new_parameters,
                                                 const Index* new_outputs_index) const=0;
  
  /**
   * @brief Evaluates the smcfixed part of likelihoods.
   *
   * @param index The index.
   */
  virtual void evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  /*
  virtual void evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                     const Parameters &conditioned_on_parameters)=0;
  */
  virtual double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  /*
  virtual double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                      const Parameters &conditioned_on_parameters)=0;
  */
  virtual double evaluate_likelihoods(const Index* index) const=0;

  // Sum of the values most recently computed for the factors in index, read
  // from each factor's cache without evaluating anything. Valid only if those
  // factors were last evaluated at the particle's current parameters.
  virtual double get_cached_likelihoods(const Index* index) const=0;

  // Consistency check for MCMC factor indices: for every factor NOT in
  // new_outputs_index whose output is exact (deterministic), re-evaluate a copy
  // at the particle's current parameters and compare with the cached value.
  // Returns the 0-based numbers of the factors whose value changed. Random
  // factors cannot be checked this way and are skipped.
  virtual std::vector<size_t> find_changed_exact_factors(const Index* new_outputs_index) const=0;
  /*
  virtual double evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed likelihoods operation.
   *
   * @param index The index.
   *
   * @return The result.
   */
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  
  /*
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)=0;
  */
  virtual double subsample_evaluate_likelihoods(const Index* index) const=0;
  /*
  virtual double subsample_evaluate_likelihoods(const Index* index,
                                                const Parameters &conditioned_on_parameters)=0;
  */
  
  /*
  virtual arma::mat direct_get_gradient_of_log(const std::string &variable)=0;
  virtual arma::mat direct_get_gradient_of_log(const std::string &variable,
                                               const Parameters &conditioned_on_parameters)=0;
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const std::string &variable)=0;
  virtual arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable) const=0;
  /*
  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable,
                                               const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable) const=0;
  
  /*
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual const Factors* get_factors() const=0;
  
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

protected:
  
  // not stored here
  /** @brief The particle. */
  Particle* particle;
  
  // not stored here

  /**
   * @brief Copies the state of another Particle into this object.
   *
   * @param another The Particle instance to copy from.
   */
  void make_copy(const FactorVariables &another);

};
}

#endif
