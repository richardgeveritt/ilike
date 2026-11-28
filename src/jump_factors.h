#ifndef JUMPFACTORS_H
#define JUMPFACTORS_H

#include <RcppArmadillo.h>
using namespace Rcpp;

#include <memory>
#include <string>
#include <vector>
#include "likelihood_factor.h"
#include "involution_step.h"

namespace ilike
{

class LikelihoodEstimator;
class LikelihoodEstimatorOutput;

// The factors of annealed reversible jump's bridges (dev/plan.md 3.22,
// "Every annealed path on the general engine", stage B). After the jump the
// state is x = (k', z') with the reverse move's auxiliary variables u'; on
// that space the start of the path is
//   f0(z', u') = pi(k, z) q(u | z) / |J(z, u)|,  (z, u) = phi(z', u'),
// (the map is its own inverse) and the end is f1(z', u') = pi(k', z') q(u' | z').
// The bridges are f0^(1-t) f1^t, as factors of the path's SMC: the move's own
// factors with power t (the model at k'), q(u' | z') with power t
// (ReverseAuxiliaryFactor) and f0 with power 1 - t (MappedBackFactor), as
// SMC's annealing from an importance proposal adds the proposal with power
// 1 - t. The path's particle carries the reverse move's auxiliary variables,
// and the indicator's value k as <k>_proposed, under names prefixed
// ilike_jump_ (so they are not taken for the model's variables).
extern const char* jump_prefix;

class ReverseAuxiliaryFactor : public LikelihoodFactor
{
public:
  ReverseAuxiliaryFactor(InvolutionFunctions* functions_in, Data* data_in);
  ReverseAuxiliaryFactor(const ReverseAuxiliaryFactor &another);
  virtual ~ReverseAuxiliaryFactor();
  Factor* duplicate() const;
  LikelihoodFactor* likelihood_factor_duplicate() const;
  // log q(u' | x') at the path's particle
  double likelihood_evaluate(const Parameters &input) const;
  arma::mat likelihood_evaluate_gradient(const std::string &variable, const Parameters &input) const;
protected:
  void specific_set_data() {}
  std::unique_ptr<InvolutionFunctions> functions;
};

class MappedBackFactor : public LikelihoodFactor
{
public:
  // inner_in: the move's own factors (owned), evaluated at the model k
  MappedBackFactor(InvolutionFunctions* functions_in, const std::string &indicator_in,
                   const std::vector<std::string> &variables_in,
                   const std::vector<LikelihoodEstimator*> &inner_in, Data* data_in);
  MappedBackFactor(const MappedBackFactor &another);
  virtual ~MappedBackFactor();
  Factor* duplicate() const;
  LikelihoodFactor* likelihood_factor_duplicate() const;
  // log f0 at the path's particle
  double likelihood_evaluate(const Parameters &input) const;
  arma::mat likelihood_evaluate_gradient(const std::string &variable, const Parameters &input) const;
protected:
  void specific_set_data() {}
  std::unique_ptr<InvolutionFunctions> functions;
  std::string indicator;
  std::vector<std::string> variables;
  std::vector<std::unique_ptr<LikelihoodEstimator>> inner;
  // (made at the first evaluation; each copy of the factor has its own)
  mutable std::vector<std::unique_ptr<LikelihoodEstimatorOutput>> inner_outputs;
};

// (stage B) the bridges' factors' powers: 1 - t for f0 and t for q(u' | x'),
// t the jump's path variable (ilike_jump_t_<j>), and 0 where it is absent
// (outside the jump's path, where the factors are not evaluated)
double jump_start_power(const Parameters &parameters, const std::string &variable);
double jump_end_power(const Parameters &parameters, const std::string &variable);
std::string jump_path_variable(size_t jump);

// the path's particle's parameters without the prefixed names (the model's),
// and with them unprefixed (xi': what the reverse move's functions read)
Parameters jump_model_parameters(const Parameters &input);
Parameters jump_reverse_parameters(const Parameters &input);

}

#endif
