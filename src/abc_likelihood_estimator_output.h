#ifndef ABCLIKELIHOODESTIMATOROUTPUT_H
#define ABCLIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include "likelihood_estimator_output.h"
#include "particles.h"

namespace ilike
{

  class ABCLikelihoodEstimator;

  class ABCLikelihoodEstimatorOutput : public LikelihoodEstimatorOutput
  {

  public:
    ABCLikelihoodEstimatorOutput();

    ABCLikelihoodEstimatorOutput(ABCLikelihoodEstimator *estimator_in);

    virtual ~ABCLikelihoodEstimatorOutput();

    ABCLikelihoodEstimatorOutput(const ABCLikelihoodEstimatorOutput &another);

    void operator=(const ABCLikelihoodEstimatorOutput &another);

    LikelihoodEstimatorOutput *duplicate() const;

    void simulate();

    void simulate(const Parameters &parameters);

    void evaluate_smcfixed_part(const Parameters &parameters);

    void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

    void subsample_simulate(const Parameters &parameters);

    void subsample_evaluate_smcfixed_part(const Parameters &parameters);

    void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

    LikelihoodEstimator *get_likelihood_estimator() const;

    arma::mat get_gradient_of_log(const std::string &variable,
                                  const Parameters &x);

    arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                            const Parameters &x);

    void forget_you_were_already_written_to_file();

    void close_ofstreams();

    void print(std::ostream &os) const;

  protected:
    // Stored in sampler.
    ABCLikelihoodEstimator *estimator;

    double distance;

    double scale_constant;

    void write_to_file(const std::string &directory_name,
                       const std::string &index = "");

    size_t compact_size() const;
  int adaptive_summary_size() const;
  void adaptive_summary(double *values) const;
  double smcadaptive_from_summary(const double *values,
                                  const Parameters &parameters) const;
  bool repoint_estimator(LikelihoodEstimator* estimator_in);
  void smcadaptive_from_summaries(const double *values,
                                  size_t stride,
                                  size_t n,
                                  const Parameters &parameters,
                                  double *results) const;
    void compact_values(double *values) const;
    void write_compact_values(const std::string &directory_name,
                              const std::string &index,
                              const double *values);

    void make_copy(const ABCLikelihoodEstimatorOutput &another);
  };
}

#endif
