#ifndef LIKELIHOODESTIMATOROUTPUT_H
#define LIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>
#include <iostream>
#include "parameters.h"

namespace ilike
{
class LikelihoodEstimator;
class HMMFactorVariables;
class VectorFactorVariables;
class AnnealedLikelihoodEstimatorOutput;
//class Data;

class LikelihoodEstimatorOutput
{
  
public:
  
  LikelihoodEstimatorOutput();
  
  virtual ~LikelihoodEstimatorOutput();
  
  //virtual void simulate(const Parameters &parameters);
  
  LikelihoodEstimatorOutput(const LikelihoodEstimatorOutput &another);
  
  void operator=(const LikelihoodEstimatorOutput &another);
  virtual LikelihoodEstimatorOutput* duplicate() const=0;
  
  // Simulate all auxilliary variables.
  virtual void simulate()=0;
  
  // Simulate all auxilliary variables.
  virtual void simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  double evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  virtual void evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  virtual void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  // Simulate all auxilliary variables.
  virtual void subsample_simulate(const Parameters &parameters)=0;
  
  // Evaluate likelihood given auxiliary variables.
  double subsample_evaluate(const Parameters &parameters);
  
  // Evaluate smcfixed part of likelihood given auxiliary variables.
  virtual void subsample_evaluate_smcfixed_part(const Parameters &parameters)=0;
  
  // Evaluate adaptive part of likelihood given auxiliary variables.
  virtual void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters)=0;
  
  virtual LikelihoodEstimator* get_likelihood_estimator() const=0;
  
  void change_data();
  void change_data(std::shared_ptr<Data> new_data);
  
  // This will end up a little bit more complicated than we have at the moment, since there is an interaction with
  // calling Simulate.
  virtual arma::mat get_gradient_of_log(const std::string &variable,
                                        const Parameters &x)=0;
  virtual arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                                  const Parameters &x)=0;
  
  virtual void forget_you_were_already_written_to_file()=0;
  
  void write(const std::string &directory_name);
  
  friend std::ostream& operator<<(std::ostream& os, const LikelihoodEstimatorOutput &p);
  
  virtual void print(std::ostream &os) const;
  
  double log_likelihood;
  
  double subsample_log_likelihood;
  
  bool write_to_file_flag;
  
  virtual void close_ofstreams()=0;
  
  // virtual double estimate_log_likelihood(const List &inputs,
  //                                        const List &auxiliary_variables) const=0;
  //
  // virtual List simulate_auxiliary_variables(const List &inputs) const=0;
  //
  // virtual void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                            const std::vector<List> &all_auxiliary_variables)=0;
  
protected:
  
  friend HMMFactorVariables;
  friend VectorFactorVariables;
  friend AnnealedLikelihoodEstimatorOutput;
  
  virtual void write_to_file(const std::string &directory_name,
                             const std::string &index = "")=0;
  
  void make_copy(const LikelihoodEstimatorOutput &another);
  
};

}

#endif

