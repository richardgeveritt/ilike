#ifndef ExactLikelihoodEstimatorOutput_H
#define ExactLikelihoodEstimatorOutput_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

class ExactLikelihoodEstimator;

class ExactLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{

public:

  ExactLikelihoodEstimatorOutput();
  ExactLikelihoodEstimatorOutput(ExactLikelihoodEstimator* estimator_in);
  virtual ~ExactLikelihoodEstimatorOutput();

  ExactLikelihoodEstimatorOutput(const ExactLikelihoodEstimatorOutput &another);
  void operator=(const ExactLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void simulate();
  
  void simulate(const Parameters &parameters);
  //double evaluate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  LikelihoodEstimator* get_likelihood_estimator() const;
  
  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                const Parameters &x);

  void print(std::ostream &os) const;

protected:

  // Stored in ModelAndAlgorithm.
  ExactLikelihoodEstimator* estimator;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;

  void make_copy(const ExactLikelihoodEstimatorOutput &another);

};

#endif
