#ifndef FACTORVARIABLES_H
#define FACTORVARIABLES_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "parameters.h"

class Particle;
class Index;
class Factors;

class FactorVariables
{

public:

  FactorVariables();
  virtual ~FactorVariables();
  
  FactorVariables(Particle* particle_in);

  FactorVariables(const FactorVariables &another);

  void operator=(const FactorVariables &another);
  virtual FactorVariables* duplicate() const=0;
  
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
  virtual double evaluate_likelihoods(const Index* index)=0;
  /*
  virtual double evaluate_likelihoods(const Index* index,
                                      const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index)=0;
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index)=0;
  
  /*
  virtual void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                               const Parameters &conditioned_on_parameters)=0;
  virtual double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                                const Parameters &conditioned_on_parameters)=0;
  */
  virtual double subsample_evaluate_likelihoods(const Index* index)=0;
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
                                               const std::string &variable)=0;
  /*
  virtual arma::mat direct_get_gradient_of_log(const Index* index,
                                               const std::string &variable,
                                               const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable)=0;
  
  /*
  virtual arma::mat direct_subsample_get_gradient_of_log(const Index* index,
                                                         const std::string &variable,
                                                         const Parameters &conditioned_on_parameters)=0;
  */
  
  virtual Factors* get_factors()=0;
  
  void set_particle(Particle* particle_in);
  Particle* get_particle();
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index) const=0;
  
  virtual void close_ofstreams()=0;

protected:
  
  // not stored here
  Particle* particle;
  
  // not stored here

  void make_copy(const FactorVariables &another);

};

#endif
