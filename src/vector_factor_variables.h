#ifndef VECTORFACTORVARIABLES_H
#define VECTORFACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include "factor_variables.h"

class VectorFactors;

class LikelihoodEstimatorOutput;

class VectorFactorVariables : public FactorVariables
{

public:

  VectorFactorVariables();
  
  VectorFactorVariables(const VectorFactors* vector_factors,
                        const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in);
  
  VectorFactorVariables(VectorFactors* vector_factors,
                        const std::vector<LikelihoodEstimatorOutput*> &likelihood_estimator_outputs_in,
                        Particle* particle_in);

  virtual ~VectorFactorVariables();

  VectorFactorVariables(const VectorFactorVariables &another);

  void operator=(const VectorFactorVariables &another);
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
  
  void forget_you_were_already_written_to_file();
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index) const;
  
  void close_ofstreams();
  
protected:

  // not stored here
  const VectorFactors* vector_factors;
  
  void make_copy(const VectorFactorVariables &another);
  
  // stored here
  std::vector<LikelihoodEstimatorOutput*> likelihood_estimator_outputs;

};

#endif
