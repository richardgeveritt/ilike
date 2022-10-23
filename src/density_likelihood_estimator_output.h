#ifndef DensityLikelihoodEstimatorOutput_H
#define DensityLikelihoodEstimatorOutput_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

class DensityLikelihoodEstimator;
class DensityEstimator;

class DensityLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{

public:

  DensityLikelihoodEstimatorOutput();
  DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in);
  virtual ~DensityLikelihoodEstimatorOutput();

  DensityLikelihoodEstimatorOutput(const DensityLikelihoodEstimatorOutput &another);
  void operator=(const DensityLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void simulate();
  
  void simulate(const Parameters &conditioned_on_parameters);
  //double evaluate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  
  void subsample_simulate(const Parameters &conditioned_on_parameters);
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
  DensityLikelihoodEstimator* estimator;
  
  // Stored here.
  //DensityEstimator* density_estimator;
  
  //Particles particles;
  std::vector<Parameters> points;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;

  void make_copy(const DensityLikelihoodEstimatorOutput &another);

};

#endif
