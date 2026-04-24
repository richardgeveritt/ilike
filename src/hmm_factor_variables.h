#ifndef HMMFACTORVARIABLES_H
#define HMMFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factor_variables.h"

namespace ilike
{
  /**
   * @file hmm_factor_variables.h
   * @brief Defines the HMMFactors class.
   *
   * Provides hmm factors functionality.
   *
   * @namespace ilike
   * @class HMMFactors
   * @brief The hmm factors class.
   */


class HMMFactors;
class Index;
class LikelihoodEstimatorOutput;
//class Data;
class Factors;

class HMMFactorVariables : public FactorVariables
{
  
public:
  
  /**
   * @brief Performs the hmmfactorvariables operation.
   */
  HMMFactorVariables();
  
  HMMFactorVariables(const HMMFactors* factors_in,
                     const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in);
  
  HMMFactorVariables(const HMMFactors* factors_in,
                     const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                     Particle* particle_in);
  
  /**
   * @brief Performs the ~hmmfactorvariables operation.
   */
  virtual ~HMMFactorVariables();
  
  /**
   * @brief Performs the hmmfactorvariables operation.
   *
   * @param another The HMMFactors instance to copy from.
   */
  HMMFactorVariables(const HMMFactorVariables &another);
  
  /**
   * @brief Assignment operator for HMMFactors.
   *
   * @param another The HMMFactors instance to copy from.
   */
  void operator=(const HMMFactorVariables &another);
  /**
   * @brief Creates a deep copy of this HMMFactors object.
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
  
  /*
   arma::mat direct_get_gradient_of_log(const std::string &variable);
   arma::mat direct_get_gradient_of_log(const std::string &variable,
   const Parameters &conditioned_on_parameters);
   
   arma::mat direct_subsample_get_gradient_of_log(const std::string &variable);
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
  
  // Stored in likelihood_estimator.
  /** @brief The hmm factors. */
  const HMMFactors* hmm_factors;
  
  /** @brief The dynamic smcfixed part. */
  double dynamic_smcfixed_part;
  
  /**
   * @brief Copies the state of another HMMFactors into this object.
   *
   * @param another The HMMFactors instance to copy from.
   */
  void make_copy(const HMMFactorVariables &another);
  
  // stored here
  //LikelihoodEstimatorOutput* initial_prior;
  
  // stored here
  // these are the likelihood terms for each measurement
  /** @brief The likelihood estimator outputs. */
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  /** @brief The likelihood estimator temp data. */
  std::vector<Data*> likelihood_estimator_temp_data;
  
};
}

#endif
