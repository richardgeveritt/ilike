#ifndef BUILTIN_REGISTRY_H
#define BUILTIN_REGISTRY_H

// The registry of built-ins (dev/plan.md 3.14, I2): one table entry per
// built-in (ilike::norm, ilike::mvnorm_rw, ilike::scaled_am, ...) giving its
// arguments (names, types, defaults) and a constructor for each role it can
// play. One function binds the arguments of a call to its entry, so every
// built-in gets the same checks, keyword arguments and named run-time
// constants, whichever front end produced the call.

#include <RcppArmadillo.h>

#include <functional>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "recipe_value.h"

namespace ilike
{

class DistributionFactor;
class IndependentProposalKernel;
class ProposalKernel;
class SymmetricProposalKernel;
class MCMCAdaptor;
class SMCAdaptor;
class MCMCTermination;
class SMCCriterion;
class SMCTermination;
class EnsembleShifter;

// What a built-in is asked to be: one value per block that takes a built-in.
enum class BuiltinRole
{
  Prior,                // prior: a DistributionFactor
  ISProposal,           // is_proposal: an IndependentProposalKernel
  IMHProposal,          // imh_proposal: an IndependentProposalKernel
  MHProposal,           // mh_proposal: a ProposalKernel
  UnadjustedProposal,   // unadjusted_proposal: a ProposalKernel
  MProposal,            // m_proposal: a SymmetricProposalKernel
  TransitionModel,      // transition_model (state-space models)
  TransitionProposal,   // transition_proposal (state-space models)
  MCMCAdaptor,          // m_adaptor, mh_adaptor
  SMCAdaptor,           // smc_adaptor
  MCMCTermination,      // mcmc_termination
  AdaptiveResampling,   // adaptive_resampling: an SMCCriterion
  AdaptiveTarget,       // adaptive_target: an SMCCriterion
  SMCTermination,       // smc_termination
  EnKShifter,           // enk_shifter
  SMCSequence,          // smc_sequence: bound only (get_smc_sequencer_info
                        // collects the sequences; the schedules are values)
  EnKSequence,          // enk_sequence (EnKI): bound only, as smc_sequence
  // bound only: the arguments are read by the code that builds the factor or
  // filter, which has the context (random-number generator, data, the user's
  // functions, the inner model) that construction needs
  SBILikelihood,        // sbi_likelihood
  AlgorithmicLikelihood,// algorithmic_likelihood (the inner model as argument "model")
  SSM,                  // ssm
  LinearGaussianDataModel,       // linear_gaussian_data_model
  LinearGaussianTransitionModel  // linear_gaussian_transition_model
};

// the block name of a role (for messages), and what it builds ("MCMC adaptor")
std::string role_block(BuiltinRole role);
std::string role_description(BuiltinRole role);

enum class ArgumentType { Variable, Real, Integer, Logical, String, Vector, Matrix, Flag };

std::string argument_type_name(ArgumentType type);

struct ArgumentSpec
{
  std::string name;
  ArgumentType type;
  bool required = true;
  // the value used when an optional argument is not given (Null: none, so
  // the constructor decides, e.g. a target acceptance rate that depends on
  // the proposal; given() then says whether it was supplied)
  RecipeValue default_value;
  std::string description;
  // bounds, checked for Real and Integer arguments
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  bool lower_open = false;
  bool upper_open = false;

  ArgumentSpec &between(double lower_in, double upper_in, bool open = true);
  ArgumentSpec &at_least(double lower_in);
  ArgumentSpec &above(double lower_in);
};

// a required argument, an optional one with a default, an optional one
// without (the constructor decides), a flag (FALSE unless given)
ArgumentSpec required(const std::string &name, ArgumentType type, const std::string &description);
ArgumentSpec optional(const std::string &name, ArgumentType type, const RecipeValue &default_value,
                      const std::string &description);
ArgumentSpec optional(const std::string &name, ArgumentType type, const std::string &description);
ArgumentSpec flag(const std::string &name, const std::string &description);

// The arguments of one call, bound to the entry's argument list: checked,
// converted to their types, defaults filled in.
class BuiltinArguments
{

public:

  BuiltinArguments(const std::string &builtin_in, BuiltinRole role_in);

  const std::string &builtin() const { return this->builtin_name; }
  BuiltinRole role() const { return this->builtin_role; }
  // "ilike::mvnorm_rw (mh_proposal)", for messages
  std::string where() const;

  bool given(const std::string &name) const;  // supplied in the call
  bool has(const std::string &name) const;    // supplied or defaulted

  std::string variable(const std::string &name) const;
  std::string string(const std::string &name) const;
  double real(const std::string &name) const;
  long integer(const std::string &name) const;
  bool logical(const std::string &name) const;
  arma::colvec vector(const std::string &name) const;
  arma::mat matrix(const std::string &name) const;

  void set(const std::string &name, const RecipeValue &value, bool given_in);

private:

  std::string builtin_name;
  BuiltinRole builtin_role;
  std::map<std::string, RecipeValue> values;
  std::set<std::string> given_names;

  const RecipeValue &get(const std::string &name) const;

};

// what an MCMC adaptor needs to know about the proposal it adapts
struct MCMCAdaptorContext
{
  std::string block;          // m_adaptor or mh_adaptor
  std::string proposal_type;  // the built-in proposal's name, "" if written by the user
};

struct Builtin
{
  std::string name;
  std::string description;
  std::vector<ArgumentSpec> arguments;
  std::vector<BuiltinRole> roles;

  // one constructor per kind of object; a role uses the one for its kind
  // (Prior: prior, and independent_proposal when the prior is used as a
  // proposal; ISProposal and IMHProposal: independent_proposal; MHProposal,
  // UnadjustedProposal: proposal, or else symmetric_proposal; MProposal,
  // TransitionModel, TransitionProposal: symmetric_proposal)
  std::function<DistributionFactor*(const BuiltinArguments&)> prior;
  std::function<IndependentProposalKernel*(const BuiltinArguments&)> independent_proposal;
  std::function<ProposalKernel*(const BuiltinArguments&)> proposal;
  std::function<SymmetricProposalKernel*(const BuiltinArguments&)> symmetric_proposal;
  std::function<MCMCAdaptor*(const BuiltinArguments&, const MCMCAdaptorContext&)> mcmc_adaptor;
  std::function<SMCAdaptor*(const BuiltinArguments&, ProposalKernel*)> smc_adaptor;
  std::function<MCMCTermination*(const BuiltinArguments&)> mcmc_termination;
  std::function<SMCCriterion*(const BuiltinArguments&, size_t number_of_particles)> smc_criterion;
  std::function<SMCTermination*(const BuiltinArguments&)> smc_termination;
  std::function<EnsembleShifter*(const BuiltinArguments&)> enk_shifter;

  bool has_role(BuiltinRole role) const;
};

// A call to a built-in in canonical form: {builtin: name, arguments: list},
// each argument named (keyword) or not (positional). Today's recipe nodes
// ({type, variables, parameters}, or {method, values} for settings passed
// from R) are rewritten into it; an argument given as text ("p1", "2",
// "scale=2") is split into its name and value here.
RecipeValue canonical_builtin_call(const RecipeValue &node);

// The entry for this name and role (an error naming the alternatives if
// there is none), and the bound arguments of a call. constants: the
// run-time constants (a list; by name, or as p1, p2, ... in old files).
const Builtin &find_builtin(const std::string &name, BuiltinRole role);
BuiltinArguments bind_builtin(const RecipeValue &node, BuiltinRole role,
                              const RecipeValue &constants);

// Construction: find, bind, and call the constructor for the role.
DistributionFactor *make_builtin_prior(const RecipeValue &node, const RecipeValue &constants);
IndependentProposalKernel *make_builtin_independent_proposal(const RecipeValue &node, BuiltinRole role,
                                                             const RecipeValue &constants);
ProposalKernel *make_builtin_proposal(const RecipeValue &node, BuiltinRole role,
                                      const RecipeValue &constants);
SymmetricProposalKernel *make_builtin_symmetric_proposal(const RecipeValue &node, BuiltinRole role,
                                                         const RecipeValue &constants);
MCMCAdaptor *make_builtin_mcmc_adaptor(const RecipeValue &node, const MCMCAdaptorContext &context,
                                       const RecipeValue &constants);
SMCAdaptor *make_builtin_smc_adaptor(const RecipeValue &node, ProposalKernel *proposal,
                                     const RecipeValue &constants);
MCMCTermination *make_builtin_mcmc_termination(const RecipeValue &node, const RecipeValue &constants);
// role: AdaptiveResampling or AdaptiveTarget
SMCCriterion *make_builtin_smc_criterion(const RecipeValue &node, BuiltinRole role,
                                         const RecipeValue &constants, size_t number_of_particles);
SMCTermination *make_builtin_smc_termination(const RecipeValue &node, const RecipeValue &constants);
EnsembleShifter *make_builtin_enk_shifter(const RecipeValue &node, const RecipeValue &constants);
// a shifter named by a string argument of another built-in (ienki_abc, enkf)
EnsembleShifter *make_builtin_enk_shifter(const std::string &name);

// the whole table (builtins.cpp)
const std::vector<Builtin> &builtin_registry();

// the table as a recipe value: one element per entry, with its roles and
// arguments (name, type, required, default, description), for the
// documentation and tests
RecipeValue builtin_table();

}

#endif
