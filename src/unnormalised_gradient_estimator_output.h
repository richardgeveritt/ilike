#ifndef UNNORMALISEDGRADIENTESTIMATOROUTPUT_H
#define UNNORMALISEDGRADIENTESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"

namespace ilike
{

class UnnormalisedGradientEstimator;

// The output of UnnormalisedGradientEstimator: value 0; the gradient estimate
// made from fields drawn at this particle's parameters, at the first request.
class UnnormalisedGradientEstimatorOutput : public LikelihoodEstimatorOutput
{
public:

  UnnormalisedGradientEstimatorOutput();
  UnnormalisedGradientEstimatorOutput(UnnormalisedGradientEstimator* estimator_in);
  virtual ~UnnormalisedGradientEstimatorOutput();

  UnnormalisedGradientEstimatorOutput(const UnnormalisedGradientEstimatorOutput &another);
  void operator=(const UnnormalisedGradientEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void simulate();
  void simulate(const Parameters &parameters);
  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  void subsample_simulate(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

  LikelihoodEstimator* get_likelihood_estimator() const;
  arma::mat get_gradient_of_log(const std::string &variable, const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable, const Parameters &x);

  void forget_you_were_already_written_to_file();
  void close_ofstreams();
  void write_to_file(const std::string &directory_name, const std::string &index = "");

protected:

  // Stored in ModelAndAlgorithm.
  UnnormalisedGradientEstimator* estimator;

  bool drawn;
  double gradient;

  void make_copy(const UnnormalisedGradientEstimatorOutput &another);
};

}

#endif
