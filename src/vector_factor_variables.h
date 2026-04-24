#ifndef VECTORFACTORVARIABLES_H
#define VECTORFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factor_variables.h"

namespace ilike
{
  /**
   * @file vector_factor_variables.h
   * @brief Defines the VectorFactors class.
   *
   * Provides vector factors functionality.
   *
   * @namespace ilike
   * @class VectorFactors
   * @brief The vector factors class.
   */


class VectorFactors;

class LikelihoodEstimatorOutput;

class VectorFactorVariables : public FactorVariables
{
  
public:
  
  /**
   * @brief Performs the vectorfactorvariables operation.
   */
  VectorFactorVariables();
  
  VectorFactorVariables(const VectorFactors* vector_factors,
                        const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in);
  
  VectorFactorVariables(VectorFactors* vector_factors,
                        const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                        Particle* particle_in);
  
  /**
   * @brief Performs the ~vectorfactorvariables operation.
   */
  virtual ~VectorFactorVariables();
  
  /**
   * @brief Performs the vectorfactorvariables operation.
   *
   * @param another The VectorFactors instance to copy from.
   */
  VectorFactorVariables(const VectorFactorVariables &another);
  
  /**
   * @brief Assignment operator for VectorFactors.
   *
   * @param another The VectorFactors instance to copy from.
   */
  void operator=(const VectorFactorVariables &another);
  /**
   * @brief Creates a deep copy of this VectorFactors object.
   *
   * @return The result.
   */
  FactorVariables* duplicate() const;
  
  /**
   * @brief Evaluates the smcfixed part of likelihoods.
   *
   * @param index The index.
   */
  void evaluate_smcfixed_part_of_likelihoods(const Index* index);
  /*
   void evaluate_smcfixed_part_of_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  /*
   double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  
  double evaluate_likelihoods(const Index* index) const;
  /*
   double evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index);
  /**
   * @brief Performs the subsample evaluate smcadaptive part given smcfixed likelihoods operation.
   *
   * @param index The index.
   *
   * @return The result.
   */
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  
  /*
   void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  
  double subsample_evaluate_likelihoods(const Index* index) const;
  /*
   double subsample_evaluate_likelihoods(const Index* index,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat direct_get_gradient_of_log(const std::string &variable) const;
  /*
   arma::mat direct_get_gradient_of_log(const std::string &variable,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable) const;
  /*
   arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat direct_get_gradient_of_log(const Index* index,
                                       const std::string &variable) const;
  /*
   arma::mat direct_get_gradient_of_log(const Index* index,
   const std::string &variable,
   const Parameters &conditioned_on_parameters);
   */
  
  arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                 const std::string &variable) const;
  /*
   arma::mat direct_subsample_get_gradient_of_log(const Index* index,
   const std::string &variable,
   const Parameters &conditioned_on_parameters);
   */
  
  const Factors* get_factors() const;
  
  /**
   * @brief Performs the forget you were already written to file operation.
   */
  void forget_you_were_already_written_to_file();
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  // not stored here
  /** @brief The vector factors. */
  const VectorFactors* vector_factors;
  
  /**
   * @brief Copies the state of another VectorFactors into this object.
   *
   * @param another The VectorFactors instance to copy from.
   */
  void make_copy(const VectorFactorVariables &another);
  
  // stored here
  /** @brief The likelihood estimator outputs. */
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;
  
};
}

#endif
