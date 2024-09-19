#ifndef ANNEALEDLIKELIHOODESTIMATOR_H
#define ANNEALEDLIKELIHOODESTIMATOR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <vector>

#include "likelihood_estimator.h"
#include "ilike_header.h"
#include "parameters.h"

namespace ilike
{
class AnnealedLikelihoodEstimatorOutput;
class IndependentProposalKernel;

/*
 class Power
 {
 
 public:
 
 Power();
 ~Power();
 
 Power(const Power &another);
 
 void operator=(const Power &another);
 Power* duplicate() const;
 
 virtual double get_power() const=0;
 
 
 
 protected:
 
 void make_copy(const Power &another);
 };
 
 class FunctionPower : public Power
 {
 
 public:
 
 FunctionPower();
 ~FunctionPower();
 
 FunctionPower(const FunctionPower &another);
 
 void operator=(const FunctionPower &another);
 FunctionPower* duplicate() const;
 
 protected:
 
 void make_copy(const FunctionPower &another);
 };
 
 class ConstantPower : public Power
 {
 
 public:
 
 DoublePower();
 ~DoublePower();
 
 DoublePower(const DoublePower &another);
 
 void operator=(const DoublePower &another);
 DoublePower* duplicate() const;
 
 protected:
 
 void make_copy(const DoublePower &another);
 };
 */


class AnnealedLikelihoodEstimator : public LikelihoodEstimator
{
  
public:
  
  AnnealedLikelihoodEstimator();
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              PowerFunctionPtr function_power_in,
                              const std::string &power_variable_in,
                              bool smcfixed_flag_in);
  
  AnnealedLikelihoodEstimator(RandomNumberGenerator* rng_in,
                              size_t* seed_in,
                              Data* data_in,
                              LikelihoodEstimator* estimator_in,
                              double constant_power_in,
                              bool smcfixed_flag_in);
  
  virtual ~AnnealedLikelihoodEstimator();
  
  AnnealedLikelihoodEstimator(const AnnealedLikelihoodEstimator &another);
  
  void operator=(const AnnealedLikelihoodEstimator &another);
  LikelihoodEstimator* duplicate() const;
  
  // double estimate_log_likelihood(const List &inputs,
  //                                const List &auxiliary_variables) const;
  
  LikelihoodEstimatorOutput* initialise();
  LikelihoodEstimatorOutput* initialise(const Parameters &parameters);
  
  void setup();
  void setup(const Parameters &parameters);
  
  // void is_setup_likelihood_estimator(const std::vector<List> &all_points,
  //                                    const std::vector<List> &all_auxiliary_variables);
  
protected:
  
  void specific_change_data(Data* new_data);
  
  friend AnnealedLikelihoodEstimatorOutput;
  
  PowerFunctionPtr function_power;
  std::string power_variable;
  double constant_power;
  bool use_constant;
  
  // Stored here.
  LikelihoodEstimator* estimator;
  
  std::ofstream log_likelihood_file_stream;
  
  // Stored here.
  //AnnealedLikelihoodEstimatorOutput* output;
  
  void make_copy(const AnnealedLikelihoodEstimator &another);
  
};
}

#endif
