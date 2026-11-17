#ifndef DensityLikelihoodEstimatorOutput_H
#define DensityLikelihoodEstimatorOutput_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;

#include <string>

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{

class DensityLikelihoodEstimator;
class DensityEstimator;
class DensityEstimatorOutput;

class DensityLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
{
  
public:
  
  DensityLikelihoodEstimatorOutput();
  DensityLikelihoodEstimatorOutput(DensityLikelihoodEstimator* estimator_in,
                                   DensityEstimator* density_estimator_in,
                                   DensityEstimator* subsample_density_estimator_in,
                                   bool store_output_in);
  virtual ~DensityLikelihoodEstimatorOutput();
  
  DensityLikelihoodEstimatorOutput(const DensityLikelihoodEstimatorOutput &another);
  void operator=(const DensityLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;
  
  void simulate();
  
  void simulate(const Parameters &conditioned_on_parameters);
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
  std::unique_ptr<DensityEstimatorOutput> density_estimator_output;
  std::unique_ptr<DensityEstimatorOutput> subsample_density_estimator_output;
  // Stored here.
  //DensityEstimator* density_estimator;
  
  std::vector<Parameters> points;
  
  double time;
  
  double log_likelihood_smcfixed_part;
  double subsample_log_likelihood_smcfixed_part;
  
  bool store_output;
  
  void write_to_file(const std::string &directory_name,
                     const std::string &index = "");

  size_t compact_size() const;
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;
  
  void compact_values(double *values) const;
  
  void write_compact_values(const std::string &directory_name,
  
                            const std::string &index,
  
                            const double *values);

  std::chrono::high_resolution_clock::time_point start_time;
  
  void make_copy(const DensityLikelihoodEstimatorOutput &another);
  
};
}

#endif
