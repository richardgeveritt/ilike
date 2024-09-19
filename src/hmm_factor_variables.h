#ifndef HMMFACTORVARIABLES_H
#define HMMFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factor_variables.h"

namespace ilike
{
class HMMFactors;
class Index;
class LikelihoodEstimatorOutput;
//class Data;
class Factors;

class HMMFactorVariables : public FactorVariables
{
  
public:
  
  HMMFactorVariables();
  
  HMMFactorVariables(const HMMFactors* factors_in,
                     const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in);
  
  HMMFactorVariables(const HMMFactors* factors_in,
                     const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                     Particle* particle_in);
  
  virtual ~HMMFactorVariables();
  
  HMMFactorVariables(const HMMFactorVariables &another);
  
  void operator=(const HMMFactorVariables &another);
  FactorVariables* duplicate() const;
  
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
  
  void forget_you_were_already_written_to_file();
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;
  
  void close_ofstreams();
  
protected:
  
  // Stored in likelihood_estimator.
  const HMMFactors* hmm_factors;
  
  double dynamic_smcfixed_part;
  
  void make_copy(const HMMFactorVariables &another);
  
  // stored here
  //LikelihoodEstimatorOutput* initial_prior;
  
  // stored here
  // these are the likelihood terms for each measurement
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;
  
  // stored here
  // data temporarily used in a likelihood estimator
  // set up to be a vector of Data* - to allow one for each llhd_estimator, but not using this funtionality at the moment - will always be one element - the same for all llhd_estimators
  std::vector<Data*> likelihood_estimator_temp_data;
  
};
}

#endif
