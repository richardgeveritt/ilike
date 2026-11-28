#ifndef CANCELLATION_PLAN_H
#define CANCELLATION_PLAN_H

// Which terms of a target cancel in an algorithm's weights and acceptance
// ratios, decided once, when the algorithm is built (dev/plan.md 3.20, S1).
//
// A weight or ratio is a product over the target's factors, at the new state
// and schedule, over the same at the old, times the backward over the forward
// kernel. A factor's terms drop out when (i) nothing it reads changes and its
// schedule does not change, or (ii) the kernel simulates exactly from it (so it
// is both a target term and the proposal density). The plan records the
// consequences for the objects the engine builds: which factors enter the
// weights, which are annealed, and whether a proposal's density is evaluated.
// It is built once per algorithm, so it costs nothing per particle.

#include <string>
#include <vector>

namespace ilike
{

// a factor of the recipe, as the target sees it
enum class FactorKind
{
  prior,      // a prior (evaluated log density of parameters)
  likelihood, // any other factor (exact, simulation-based, algorithmic, ...)
  none        // a factor entry that makes no estimator (a prior's simulator alone)
};

// the target of an algorithm: its factors, and its SMC sequence if any
struct TargetSpecification
{
  std::vector<FactorKind> kinds;  // per recipe factor
  std::vector<bool> listed;       // per recipe factor: named by the sequence
  bool sequence = false;          // the algorithm has an SMC sequence
  bool annealing = false;         // ... and it is annealing (tempering)
  std::string annealing_variable; // the sequence variable that is the power
  bool factors_listed = false;    // the sequence names the factors it affects
  bool data_tempering = false;    // the sequence is data tempering (S2)
  std::string data_variable;      // ... and its variable (the fraction of the data)
  size_t number_of_sequences = 0;
};

// the options of a data tempering sequence, data(variable, schedule, factors,
// order, anneal) (dev/plan.md 3.20, S2)
struct DataTemperingOptions
{
  bool random_order = true;       // order = random (else: the order given)
  bool anneal = true;             // anneal each data point in (else: add whole data points)
};

// what a proposal or kernel simulates exactly from (rule ii)
struct KernelDeclaration
{
  bool simulates_priors = false;           // an importance proposal that is the prior
  bool simulates_transition_model = false; // a transition proposal that is the transition model
};

// the importance proposal is the prior
inline KernelDeclaration prior_as_proposal()
{
  KernelDeclaration declaration;
  declaration.simulates_priors = true;
  return declaration;
}

// the transition proposal is the transition model (the bootstrap filter)
inline KernelDeclaration transition_model_as_proposal()
{
  KernelDeclaration declaration;
  declaration.simulates_transition_model = true;
  return declaration;
}

// what happens to one factor
struct FactorPlan
{
  bool in_weight = false;         // enters importance and SMC weights
  bool in_without_priors = false; // enters the index of likelihood factors
  bool annealed = false;          // raised to the annealing power
  bool fixed_in_smc = true;       // its value does not change along the sequence (only a power may)
  bool data_tempered = false;     // its data points enter along the sequence (S2)
};

struct CancellationPlan
{
  std::vector<FactorPlan> factors;       // per recipe factor
  bool annealing = false;
  std::string annealing_variable;
  bool add_proposal_path_factor = false; // the q^(1 - p) factor that starts an annealing path at q
  // (2026-11-27) the positions of the factors the algorithm appends after the
  // model's (q^(1 - p); annealed reversible jump's bridges): they read the
  // moved variables, so they join every move's factor index
  std::vector<size_t> appended_factors;
  std::string data_variable;             // (S2) the variable of data tempering
  bool proposal_is_evaluated = true;     // the importance proposal's density enters the first weights
  std::vector<size_t> followed_observations; // (3.21 (c2'')) the observations of each filter following data tempering
};

// a kernel's density enters a weight unless the kernel simulates exactly from
// factors of the target, which then cancel with it (rule ii)
inline bool kernel_density_is_evaluated(const KernelDeclaration &kernel)
{
  return !(kernel.simulates_priors || kernel.simulates_transition_model);
}

// the plan for a target whose importance proposal is declared by
// importance_proposal (a default declaration when there is no importance
// proposal, as in MCMC); stops with a message for a prior listed for annealing
// that cancels
CancellationPlan plan_cancellations(const TargetSpecification &target,
                                    const KernelDeclaration &importance_proposal);

}

#endif
