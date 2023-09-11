#ifndef DensityLikelihoodEstimatorOutput_H
#define DensityLikelihoodEstimatorOutput_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <string>

#include "likelihood_estimator_output.h"
#include "particles.h"

class DensityLikelihoodEstimator;
class DensityEstimator;
class DensityEstimatorOutput;

class DensityLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{

public:

  DensityLikelihoodEstimatorOutput();
  DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
                                   DensityEstimator* density_estimator_in,
                                   DensityEstimator* subsample_density_estimator_in);
  //DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
  //                                 const Parameters &conditioned_on_parameters);
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
  
  void fit(const std::vector<Parameters> &points);
  void subsample_fit(const std::vector<Parameters> &points);

  void forget_you_were_already_written_to_file();
  
  void close_ofstreams();
  
  void print(std::ostream &os) const;

protected:

  // Stored in ModelAndAlgorithm.
  DensityLikelihoodEstimator* estimator;
  
  // Stored here.
  DensityEstimatorOutput* density_estimator_output;
  DensityEstimatorOutput* subsample_density_estimator_output;
  // Stored here.
  //DensityEstimator* density_estimator;
  
  //Particles particles;
  std::vector<Parameters> points;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");
  
  //std::string results_name;

  void make_copy(const DensityLikelihoodEstimatorOutput &another);

};

#endif
