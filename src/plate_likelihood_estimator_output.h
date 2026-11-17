#ifndef PLATELIKELIHOODESTIMATOROUTPUT_H
#define PLATELIKELIHOODESTIMATOROUTPUT_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <cstdint>
#include <vector>
#include "exact_likelihood_estimator_output.h"

namespace ilike
{

class PlateLikelihoodEstimator;

// The output of a plate factor (dev/plan.md 3.6 Layer 1, I4b): the value of
// each term at this particle's parameters. When a move makes the output of a
// proposed particle, carry_over_from copies the current particle's terms and
// marks the terms that read an element whose value differs between the two
// parameter sets; the next evaluation recomputes only those. Every other way
// of reaching a state (reinitialise, restoring from a summary or a compact
// record) leaves no valid terms, so the next evaluation recomputes them all.
// The factor's value is always the sum of all the terms in order, so a
// partial evaluation gives exactly the value of a full one.
class PlateLikelihoodEstimatorOutput : public ExactLikelihoodEstimatorOutput
{
public:

  PlateLikelihoodEstimatorOutput();
  PlateLikelihoodEstimatorOutput(PlateLikelihoodEstimator* estimator_in);
  virtual ~PlateLikelihoodEstimatorOutput();

  PlateLikelihoodEstimatorOutput(const PlateLikelihoodEstimatorOutput &another);
  void operator=(const PlateLikelihoodEstimatorOutput &another);
  LikelihoodEstimatorOutput* duplicate() const;

  void evaluate_smcfixed_part(const Parameters &parameters);
  void evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);
  void subsample_evaluate_smcfixed_part(const Parameters &parameters);
  void subsample_evaluate_smcadaptive_part_given_smcfixed(const Parameters &parameters);

  arma::mat get_gradient_of_log(const std::string &variable,
                                const Parameters &x);
  arma::mat subsample_get_gradient_of_log(const std::string &variable,
                                          const Parameters &x);

  void carry_over_from(const LikelihoodEstimatorOutput &previous,
                       const Parameters &previous_parameters,
                       const Parameters &new_parameters);

protected:

  bool reinitialise(const Parameters &parameters);
  bool assign_from(const LikelihoodEstimatorOutput &another);
  bool repoint_estimator(LikelihoodEstimator* estimator_in);

  double plate_evaluate(const Parameters &parameters);

  // Stored in ModelAndAlgorithm.
  PlateLikelihoodEstimator* plate;

  std::vector<double> term_values;
  // term_values hold the terms at this output's parameters (except those
  // listed in pending, when partial), for the plate's structure_version
  bool terms_valid;
  size_t terms_version;
  bool partial;
  std::vector<uint32_t> pending;
  std::vector<char> marked;

  void make_copy(const PlateLikelihoodEstimatorOutput &another);
};

}

#endif
