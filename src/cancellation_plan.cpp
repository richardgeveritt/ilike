#include <Rcpp.h>
#include "cancellation_plan.h"

namespace ilike
{

CancellationPlan plan_cancellations(const TargetSpecification &target,
                                    const KernelDeclaration &importance_proposal)
{
  CancellationPlan plan;
  plan.annealing = target.annealing;
  plan.annealing_variable = target.annealing_variable;
  plan.data_variable = target.data_variable;
  plan.factors.resize(target.kinds.size());

  // (S2) data tempering acts alone, on likelihood factors
  if (target.data_tempering && (target.number_of_sequences > 1))
    Rcpp::stop("smc_sequence: data tempering (data(...)) cannot be combined with another sequence yet.");

  // rule (ii): an importance proposal that is the prior simulates exactly
  // from the priors, so they cancel with it in the weights
  const bool priors_cancel = importance_proposal.simulates_priors;

  for (size_t i = 0; i < target.kinds.size(); ++i)
  {
    FactorPlan &factor = plan.factors[i];
    const bool listed = target.listed[i];
    // annealing raises the listed factors to the power (all factors if none
    // are listed); the others are unchanged along the sequence, so they
    // cancel between successive targets (rule i)
    const bool annealed_if_in_target = target.annealing && (!target.factors_listed || listed);
    // another sequence (a parameter of the factors) changes the listed
    // factors (all, if none are listed): those are evaluated again when the
    // target changes; under annealing a factor's value is fixed and only its
    // power changes (the cost saving of evaluating it once)
    const bool changes_with_sequence = target.sequence && (!target.factors_listed || listed);
    factor.fixed_in_smc = target.annealing || !changes_with_sequence;

    switch (target.kinds[i])
    {
      case FactorKind::prior:
        if (target.data_tempering && target.factors_listed && listed)
          Rcpp::stop("smc_sequence data tempering: factor " + std::to_string(i + 1) +
                     " is a prior; data tempering applies to likelihood factors.");
        if (priors_cancel && target.annealing && listed)
          Rcpp::stop("smc_sequence annealing: factor " + std::to_string(i + 1) +
                     " is a prior that is also used as the importance proposal, so it "
                     "cancels and cannot be annealed. Remove it from the factors to "
                     "anneal, or give an is_proposal.");
        factor.in_weight = !priors_cancel;
        factor.in_without_priors = false;
        if (target.data_tempering)
          factor.fixed_in_smc = true; // (data tempering leaves the priors unchanged)
        factor.annealed = annealed_if_in_target && !priors_cancel;
        break;
      case FactorKind::likelihood:
        factor.in_weight = true;
        factor.in_without_priors = true;
        factor.annealed = annealed_if_in_target;
        // (S2) the listed factors (all the likelihood factors, if none are
        // listed) take their data points in along the sequence; their value
        // changes with it, so they are evaluated again at each step
        factor.data_tempered = target.data_tempering && (!target.factors_listed || listed);
        break;
      case FactorKind::none:
        break;
    }
  }

  // annealing every factor from an importance proposal q that is not the
  // prior: the path (all factors)^p q^(1 - p) starts at q, so q enters as a
  // factor with power 1 - p; when only chosen factors are annealed, the path
  // is (other factors) (chosen factors)^p and the first weights are target / q
  plan.add_proposal_path_factor = !priors_cancel && target.annealing && !target.factors_listed;

  // rule (ii) again: the prior as proposal is not evaluated in the first weights
  plan.proposal_is_evaluated = kernel_density_is_evaluated(importance_proposal);
  return plan;
}

}
