#include "recipe_internal.h"

IndependentProposalKernel *
get_proposal(const RecipeValue &model, const RecipeValue &model_parameters, Data *data) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.contains("is_proposal")) {
    RecipeValue is_proposals = model["is_proposal"];

    for (size_t i = 0; i < is_proposals.size(); ++i) {
      if ((is_proposals[i]).is_list()) {
        RecipeValue current_proposal = is_proposals[i];
        if (current_proposal.contains("is_proposal")) {
          RecipeValue proposal_info = current_proposal["is_proposal"];

          if ((proposal_info).is_list()) {
            if (proposal_info.contains("type") &&
                proposal_info.contains("variables")) {
              proposal = make_builtin_independent_proposal(proposal_info, BuiltinRole::ISProposal, model_parameters);
            } else {
              stop("Missing information for proposal in model file.");
            }
          } else {
            stop("Error in proposal section of model file.");
          }
        } else if (current_proposal.contains(
                       "evaluate_log_is_proposal") &&
                   current_proposal.contains(
                       "simulate_is_proposal") &&
                   current_proposal.contains("type")) {
          RecipeValue evaluate_log_is_proposal_value =
              current_proposal["evaluate_log_is_proposal"];
          RecipeValue simulate_is_proposal_value =
              current_proposal["simulate_is_proposal"];
          size_t type = current_proposal["type"];
          if (type == 1) {
            proposal = new CustomDistributionProposalKernel(
                (simulate_is_proposal_value).function_pointer<SimulateDistributionPtr>(),
                (evaluate_log_is_proposal_value).function_pointer<EvaluateLogDistributionPtr>());
          } else if (type == 2) {
            proposal = new CustomIndependentProposalKernel(
                (simulate_is_proposal_value).function_pointer<SimulateIndependentProposalPtr>(),
                (evaluate_log_is_proposal_value).function_pointer<EvaluateLogIndependentProposalPtr>());
          } else if (type == 3) {
            proposal = new CustomGuidedDistributionProposalKernel(
                (simulate_is_proposal_value).function_pointer<SimulateGuidedDistributionPtr>(),
                (evaluate_log_is_proposal_value).function_pointer<EvaluateLogGuidedDistributionPtr>(),
                data);
          } else if (type == 4) {
            proposal = new CustomGuidedIndependentProposalKernel(
                (simulate_is_proposal_value).function_pointer<SimulateGuidedIndependentProposalPtr>(),
                (evaluate_log_is_proposal_value).function_pointer<EvaluateLogGuidedIndependentProposalPtr>(),
                data);
          } else {
            stop("Functions read in for proposal are invalid.");
          }
        } else {
          stop("Invalid importance proposal. Maybe you specified a method to "
               "simulate from the proposal, but no method to evaluate it?");
        }
      } else {
        stop("Error in importance proposals section of model file.");
      }

      if (proposal != NULL) {
        proposals.push_back(proposal);
      }
    }
  }

  if (proposal == NULL) {
    stop("No proposal specified.");
  } else {
    if (proposals.size() == 1) {
      return proposal;
    } else {
      return new CompositeIndependentProposalKernel(proposals);
    }
  }
}

IndependentProposalKernel *get_prior_as_proposal(const RecipeValue &model,
                                                 const RecipeValue &model_parameters) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {
      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];
        if (current_factor.contains("prior")) {
          if ((current_factor["prior"]).is_list()) {
            RecipeValue current_prior = current_factor["prior"];
            if (current_prior.contains("type") &&
                current_prior.contains("variables")) {
              proposal = make_builtin_independent_proposal(current_prior, BuiltinRole::Prior, model_parameters);
              proposals.push_back(proposal);
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if ((current_factor.contains(
                       "evaluate_log_prior")) &&
                   (current_factor.contains("simulate_prior"))) {
          RecipeValue evaluate_log_prior_value = current_factor["evaluate_log_prior"];

          RecipeValue simulate_prior_value = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(
              (simulate_prior_value).function_pointer<SimulateDistributionPtr>(),
              (evaluate_log_prior_value).function_pointer<EvaluateLogDistributionPtr>());

          if (proposal != NULL) {
            proposals.push_back(proposal);
          }
        }
      } else {
        stop("Error in factors section of model file.");
      }
    }
  }

  if (proposal == NULL) {
    stop("No suitable prior to use as proposal.");
  }

  if (proposals.size() == 1) {
    return proposal;
  } else {
    return new CompositeIndependentProposalKernel(proposals);
  }
}

IndependentProposalKernel *
get_prior_as_simulate_only_proposal(const RecipeValue &model,
                                    const RecipeValue &model_parameters) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.contains("factor")) {
    RecipeValue factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {
      if ((factors[i]).is_list()) {
        RecipeValue current_factor = factors[i];
        if (current_factor.contains("prior")) {
          if ((current_factor["prior"]).is_list()) {
            RecipeValue current_prior = current_factor["prior"];
            if (current_prior.contains("type") &&
                current_prior.contains("variables")) {
              proposal = make_builtin_independent_proposal(current_prior, BuiltinRole::Prior, model_parameters);
              proposals.push_back(proposal);
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.contains("simulate_prior")) {
          RecipeValue simulate_prior_value = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(
              (simulate_prior_value).function_pointer<SimulateDistributionPtr>());

          if (proposal != NULL) {
            proposals.push_back(proposal);
          }
        }
      } else {
        stop("Error in factors section of model file.");
      }
    }
  }

  if (proposal == NULL) {
    stop("No suitable prior to use as proposal.");
  } else {
    if (proposals.size() == 1) {
      return proposal;
    } else {
      return new CompositeIndependentProposalKernel(proposals);
    }
  }
}

ProposalKernel *get_mh_proposal(const RecipeValue &current_proposal,
                                const RecipeValue &model_parameters, Data *data) {
  ProposalKernel *proposal;
  if (current_proposal.contains("mh_proposal")) {
    RecipeValue proposal_info = current_proposal["mh_proposal"];

    if ((proposal_info).is_list()) {
      if (proposal_info.contains("type") &&
          proposal_info.contains("variables")) {
        proposal = make_builtin_proposal(proposal_info, BuiltinRole::MHProposal, model_parameters);
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.contains(
                 "evaluate_log_mh_proposal") &&
             current_proposal.contains("simulate_mh_proposal") &&
             current_proposal.contains("type")) {
    RecipeValue evaluate_log_mh_proposal_value =
        current_proposal["evaluate_log_mh_proposal"];
    RecipeValue simulate_mh_proposal_value = current_proposal["simulate_mh_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsProposalKernel(
          (simulate_mh_proposal_value).function_pointer<SimulateNoParamsMCMCProposalPtr>(),
          (evaluate_log_mh_proposal_value).function_pointer<EvaluateLogNoParamsMCMCProposalPtr>());
    } else if (type == 2) {
      proposal = new CustomProposalKernel(
          (simulate_mh_proposal_value).function_pointer<SimulateMCMCProposalPtr>(),
          (evaluate_log_mh_proposal_value).function_pointer<EvaluateLogMCMCProposalPtr>());
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsProposalKernel(
          (simulate_mh_proposal_value).function_pointer<SimulateGuidedNoParamsMCMCProposalPtr>(),
          (evaluate_log_mh_proposal_value).function_pointer<EvaluateLogGuidedNoParamsMCMCProposalPtr>(),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedProposalKernel(
          (simulate_mh_proposal_value).function_pointer<SimulateGuidedMCMCProposalPtr>(),
          (evaluate_log_mh_proposal_value).function_pointer<EvaluateLogGuidedMCMCProposalPtr>(),
          data);
    } else {
      stop("Functions read in for proposal are invalid.");
    }
  } else {
    stop("Invalid mh proposal. Maybe you specified a method to simulate from "
         "the proposal, but no method to evaluate it?");
  }

  if (proposal == NULL) {
    stop("No suitable proposal found in model file.");
  } else {
    return proposal;
  }
}

ProposalKernel *get_unadjusted_proposal(const RecipeValue &current_proposal,
                                        const RecipeValue &model_parameters,
                                        Data *data) {
  ProposalKernel *proposal;
  if (current_proposal.contains("unadjusted_proposal")) {
    RecipeValue proposal_info = current_proposal["unadjusted_proposal"];

    if ((proposal_info).is_list()) {
      if (proposal_info.contains("type") &&
          proposal_info.contains("variables")) {
        proposal = make_builtin_proposal(proposal_info, BuiltinRole::UnadjustedProposal, model_parameters);
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.contains(
                 "simulate_unadjusted_proposal") &&
             current_proposal.contains("type")) {
    RecipeValue simulate_unadjusted_proposal_value =
        current_proposal["simulate_unadjusted_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsProposalKernel(
          (simulate_unadjusted_proposal_value).function_pointer<SimulateNoParamsMCMCProposalPtr>());
    } else if (type == 2) {
      proposal = new CustomProposalKernel(
          (simulate_unadjusted_proposal_value).function_pointer<SimulateMCMCProposalPtr>());
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsProposalKernel(
          (simulate_unadjusted_proposal_value).function_pointer<SimulateGuidedNoParamsMCMCProposalPtr>(),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedProposalKernel(
          (simulate_unadjusted_proposal_value).function_pointer<SimulateGuidedMCMCProposalPtr>(),
          data);
    } else {
      stop("Functions read in for proposal are invalid.");
    }
  } else {
    stop("Invalid mh proposal. Maybe you specified a method to simulate from "
         "the proposal, but no method to evaluate it?");
  }

  if (proposal == NULL) {
    stop("No suitable proposal found in model file.");
  } else {
    return proposal;
  }
}

IndependentProposalKernel *get_imh_proposal(const RecipeValue &current_proposal,
                                            const RecipeValue &model_parameters,
                                            Data *data) {
  IndependentProposalKernel *proposal;
  if (current_proposal.contains("imh_proposal")) {
    RecipeValue proposal_info = current_proposal["imh_proposal"];

    if ((proposal_info).is_list()) {
      if (proposal_info.contains("type") &&
          proposal_info.contains("variables")) {
        proposal = make_builtin_independent_proposal(proposal_info, BuiltinRole::IMHProposal, model_parameters);
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.contains(
                 "evaluate_log_imh_proposal") &&
             current_proposal.contains("simulate_imh_proposal") &&
             current_proposal.contains("type")) {
    RecipeValue evaluate_log_imh_proposal_value =
        current_proposal["evaluate_log_imh_proposal"];
    RecipeValue simulate_imh_proposal_value = current_proposal["simulate_imh_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomDistributionProposalKernel(
          (simulate_imh_proposal_value).function_pointer<SimulateDistributionPtr>(),
          (evaluate_log_imh_proposal_value).function_pointer<EvaluateLogDistributionPtr>());
    } else if (type == 2) {
      proposal = new CustomIndependentProposalKernel(
          (simulate_imh_proposal_value).function_pointer<SimulateIndependentProposalPtr>(),
          (evaluate_log_imh_proposal_value).function_pointer<EvaluateLogIndependentProposalPtr>());
    } else if (type == 3) {
      proposal = new CustomGuidedDistributionProposalKernel(
          (simulate_imh_proposal_value).function_pointer<SimulateGuidedDistributionPtr>(),
          (evaluate_log_imh_proposal_value).function_pointer<EvaluateLogGuidedDistributionPtr>(),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedIndependentProposalKernel(
          (simulate_imh_proposal_value).function_pointer<SimulateGuidedIndependentProposalPtr>(),
          (evaluate_log_imh_proposal_value).function_pointer<EvaluateLogGuidedIndependentProposalPtr>(),
          data);
    } else {
      stop("Functions read in for proposal are invalid.");
    }
  } else {
    stop("Invalid imh proposal. Maybe you specified a method to simulate from "
         "the proposal, but no method to evaluate it?");
  }

  if (proposal == NULL) {
    stop("No suitable prior to use as proposal.");
  } else {
    return proposal;
  }
}

SymmetricProposalKernel *get_m_proposal(const RecipeValue &current_proposal,
                                        const RecipeValue &model_parameters,
                                        Data *data) {
  SymmetricProposalKernel *proposal;
  if (current_proposal.contains("m_proposal")) {
    RecipeValue proposal_info = current_proposal["m_proposal"];

    if ((proposal_info).is_list()) {
      if (proposal_info.contains("type") &&
          proposal_info.contains("variables")) {
        proposal = make_builtin_symmetric_proposal(proposal_info, BuiltinRole::MProposal, model_parameters);
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.contains("simulate_m_proposal") &&
             current_proposal.contains("type")) {
    RecipeValue simulate_m_proposal_value = current_proposal["simulate_m_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsSymmetricProposalKernel(
          (simulate_m_proposal_value).function_pointer<SimulateNoParamsMCMCProposalPtr>());
    } else if (type == 2) {
      proposal = new CustomSymmetricProposalKernel(
          (simulate_m_proposal_value).function_pointer<SimulateMCMCProposalPtr>());
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsSymmetricProposalKernel(
          (simulate_m_proposal_value).function_pointer<SimulateGuidedNoParamsMCMCProposalPtr>(),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedSymmetricProposalKernel(
          (simulate_m_proposal_value).function_pointer<SimulateGuidedMCMCProposalPtr>(), data);
    } else {
      stop("Functions read in for proposal are invalid.");
    }
  } else {
    stop("Invalid m proposal. Maybe you specified a method to simulate from "
         "the proposal, but no method to evaluate it?");
  }

  if (proposal == NULL) {
    stop("No suitable proposal found in model file.");
  } else {
    return proposal;
  }
}

MCMCTermination *make_mcmc_termination(const RecipeValue &model_parameters,
                                       const RecipeValue &mcmc_termination_method) {
  if (!(mcmc_termination_method.contains("method") || mcmc_termination_method.contains("type")))
    stop("No valid method found for MCMC termination.");
  return make_builtin_mcmc_termination(mcmc_termination_method, model_parameters);
}

static MoveFactors set_move_factor_indices(MCMC *mcmc, const RecipeValue &proposal_info,
                                           const std::string &prefix,
                                           const Index *full_index,
                                           bool check_factor_index,
                                           size_t block_number) {
  MoveFactors move_factors;
  move_factors.label = prefix + "_proposal " + std::to_string(block_number);
  move_factors.unadjusted = (prefix == "unadjusted");
  if (full_index != NULL)
    move_factors.used.assign(full_index->begin(), full_index->end());

  std::string used_name = prefix + "_factor_index";
  std::string cut_name = prefix + "_cut_factor_index";
  bool has_used = proposal_info.contains(used_name.c_str());
  bool has_cut = proposal_info.contains(cut_name.c_str());
  if (!has_used && !has_cut)
    return move_factors;
  mcmc->set_check_factor_index(check_factor_index, move_factors.label);
  if (full_index == NULL)
    Rcpp::stop("make_mcmc: full index of factors not available.");

  std::vector<size_t> all_factors(full_index->begin(), full_index->end());

  auto read_index = [&](const std::string &name) {
    RecipeValue index_value = proposal_info[name];
    std::vector<double> index = Rcpp::as<std::vector<double>>((index_value).function_pointer<NumericVectorPtr>()());
    std::vector<size_t> result;
    for (size_t k = 0; k < index.size(); ++k) {
      double value = index[k];
      if ((value < 1) || (value != std::floor(value)) ||
          (std::find(all_factors.begin(), all_factors.end(),
                     size_t(value) - 1) == all_factors.end()))
        Rcpp::stop(name + ": " + std::to_string(value) +
                   " is not the number of a factor in the model (factors are "
                   "numbered from 1 in the order they appear).");
      if (std::find(result.begin(), result.end(), size_t(value) - 1) ==
          result.end())
        result.push_back(size_t(value) - 1);
    }
    return result;
  };

  std::vector<size_t> cut;
  if (has_cut)
    cut = read_index(cut_name);

  std::vector<size_t> used;
  if (has_used) {
    used = read_index(used_name);
  } else {
    for (auto i = all_factors.begin(); i != all_factors.end(); ++i)
      if (std::find(cut.begin(), cut.end(), *i) == cut.end())
        used.push_back(*i);
  }

  for (auto i = cut.begin(); i != cut.end(); ++i)
    if (std::find(used.begin(), used.end(), *i) != used.end())
      Rcpp::stop("Factor " + std::to_string(*i + 1) + " is in both " +
                 used_name + " and " + cut_name + ".");

  VectorIndex used_index(used); // the setters copy their arguments
  mcmc->set_index_if_null(&used_index);
  if (has_cut) {
    VectorIndex cut_index(cut);
    mcmc->set_cut_index(&cut_index);
  }
  move_factors.used = used;
  move_factors.cut = cut;
  return move_factors;
}

// Unadjusted moves do not evaluate any factor (often the reason for using them
// is to avoid this, and the user may not even supply the functions to do it).
// So after an unadjusted move, the cached values of the factors it changes
// (its factor index and cut factor index) are not valid. A Metropolis or
// Metropolis-Hastings move reads the cached values of the factors in its own
// factor index for its acceptance-ratio denominator, so it may not use any
// factor that an unadjusted move changes.
static void check_unadjusted_compatibility(const std::vector<MoveFactors> &moves) {
  for (auto u = moves.begin(); u != moves.end(); ++u) {
    if (!u->unadjusted)
      continue;
    std::vector<size_t> changed = u->used;
    changed.insert(changed.end(), u->cut.begin(), u->cut.end());
    for (auto m = moves.begin(); m != moves.end(); ++m) {
      if (m->unadjusted)
        continue;
      std::string shared;
      for (auto f = m->used.begin(); f != m->used.end(); ++f)
        if (std::find(changed.begin(), changed.end(), *f) != changed.end())
          shared = shared + (shared.empty() ? "" : ", ") + std::to_string(*f + 1);
      if (!shared.empty())
        Rcpp::stop("The unadjusted move " + u->label +
                   " does not evaluate factors, but " + m->label +
                   " uses factor(s) " + shared +
                   " (which the unadjusted move changes) in its acceptance "
                   "ratio. Unadjusted moves can only be combined with moves "
                   "whose factor index does not include the factors of the "
                   "unadjusted move (its factor index and cut factor index).");
    }
  }
}

// Adaptive MCMC (dev/plan.md 3.2): what adapts the proposal just built.
// Either a built-in adaptor block after the proposal (/***mh_adaptor,...***/,
// stored in the proposal's part of the recipe under "<prefix>_adaptor"; the
// adaptors, their arguments and defaults are in the registry, builtins.cpp),
// or the user's adaptation function (under "adapt_<prefix>_proposal"), with
// which the only block allowed is ilike::averaged() (Polyak & Juditsky 1992
// averaging of the iterates).
// The scale a proposal moves on (dev/plan.md 3.18 H2, C209): the move's
// transform = …, lowered to a list of {variable, kind, lower, upper}.
Transform *transform_from_recipe(const RecipeValue &list) {
  std::vector<ElementwiseTransformComponent::Spec> specs;
  for (size_t k = 0; k < list.size(); ++k) {
    const RecipeValue &item = list[k];
    ElementwiseTransformComponent::Spec spec;
    spec.variable = item["variable"].as_string();
    std::string kind = item["kind"].as_string();
    spec.kind = (kind == "log") ? ElementwiseTransformComponent::Kind::log
              : (kind == "lower") ? ElementwiseTransformComponent::Kind::lower
              : (kind == "upper") ? ElementwiseTransformComponent::Kind::upper
              : ElementwiseTransformComponent::Kind::logit;
    spec.lower = item["lower"].as_double();
    spec.upper = item["upper"].as_double();
    specs.push_back(spec);
  }
  if (specs.empty())
    return NULL;
  return new Transform(new ElementwiseTransformComponent(specs));
}

static void attach_transform(ProposalKernel *proposal, const RecipeValue &proposal_info) {
  if (!proposal_info.contains("transform"))
    return;
  Transform *transform = transform_from_recipe(proposal_info["transform"]);
  if (transform != NULL)
    proposal->set_transform(transform);
}

void attach_mcmc_adaptor(ProposalKernel *proposal, const RecipeValue &proposal_info,
                         const std::string &prefix,
                         const RecipeValue &model_parameters) {
  std::string adaptor_name = prefix + "_adaptor";
  std::string function_name = "adapt_" + prefix + "_proposal";
  bool has_block = proposal_info.contains(adaptor_name.c_str());
  bool has_function = proposal_info.contains(function_name.c_str());
  if (has_function) {
    // the only block allowed with a user's function: ilike::averaged()
    bool function_averaging = false;
    if (has_block) {
      RecipeValue block_info = proposal_info[adaptor_name];
      std::string block_type = block_info["type"];
      if (!iequals(block_type, "averaged"))
        stop("A proposal with an " + function_name + " function can have only "
             "the block " + adaptor_name + ",ilike::averaged() with it.");
      function_averaging = true;
    }
    RecipeValue adapt_value = proposal_info[function_name];
    MCMCAdaptor* adaptor = new ProposalParameterAdaptor(
        (adapt_value).function_pointer<AdaptProposalPtr>(), 0.6, 1.0, function_averaging);
    if (proposal_info.contains((adaptor_name + "_until").c_str()))
      adaptor->set_until(AdaptationCriterion::from_recipe(proposal_info[adaptor_name + "_until"]));
    proposal->set_mcmc_adaptor(adaptor);
    return;
  }
  if (!has_block)
    return;

  std::string proposal_type = "";
  std::string proposal_name = prefix + "_proposal";
  if (proposal_info.contains(proposal_name.c_str())) {
    RecipeValue proposal_value = proposal_info[proposal_name];
    if ((proposal_value).is_list() && proposal_value.contains("type"))
      proposal_type = (proposal_value["type"]).as_string();
  }

  RecipeValue adaptor_info = proposal_info[adaptor_name];
  MCMCAdaptorContext context;
  context.block = adaptor_name;
  context.proposal_type = proposal_type;
  MCMCAdaptor* adaptor = make_builtin_mcmc_adaptor(adaptor_info, context, model_parameters);
  // when to stop adapting (until = …, dev/plan.md 3.18 item 10)
  if (proposal_info.contains((adaptor_name + "_until").c_str()))
    adaptor->set_until(AdaptationCriterion::from_recipe(proposal_info[adaptor_name + "_until"]));
  proposal->set_mcmc_adaptor(adaptor);
}

// Adaptation of a proposal from the SMC particles (dev/plan.md 3.2, B2): the
// block /***smc_adaptor,ilike::gaussian(...)***/ after the proposal (in the
// registry, builtins.cpp), stored in the proposal's part of the recipe under "smc_adaptor". Before
// each move of an SMC algorithm (or ensemble Kalman method), the proposal's
// covariance is set to the weighted covariance of the particles (and a
// mirror's centre to their weighted mean; a principal-components proposal's
// directions to its eigenvectors); with c, the scale to c^2/d (c^2 for
// componentwise and principal-components proposals).
void attach_smc_adaptor(ProposalKernel *proposal, const RecipeValue &proposal_info,
                        const RecipeValue &model_parameters) {
  if (!proposal_info.contains("smc_adaptor"))
    return;
  proposal->set_smc_adaptor(make_builtin_smc_adaptor(proposal_info["smc_adaptor"], proposal,
                                                     model_parameters));
}

MCMC *make_mcmc(const RecipeValue &model, const RecipeValue &model_parameters, Data *data,
                const RecipeValue &mcmc_termination_method,
                const RecipeValue &mcmc_weights_method, Index *full_index) {
  MCMCTermination *termination_method =
      make_mcmc_termination(model_parameters, mcmc_termination_method);

  MCMC *mcmc = NULL;

  // Optional consistency check of the factor indices, set by the R wrappers
  // (MCMC(..., check_factor_index = TRUE)) as an element of the recipe.
  bool check_factor_index = false;
  if (model.contains("check_factor_index"))
    check_factor_index = (model["check_factor_index"]).as_bool();

  // MCMC(..., thin): store every thin-th iteration (recipe element mcmc_thin).
  size_t mcmc_thin = 1;
  if (model.contains("mcmc_thin"))
    mcmc_thin = (model["mcmc_thin"]).as_size_t();

  if (model.contains("order_of_mcmc")) {
    std::vector<double> order_of_mcmc = model["order_of_mcmc"];
    std::vector<MCMC *> moves;
    moves.reserve(order_of_mcmc.size());
    std::vector<MoveFactors> move_factors;

    size_t simulate_mh_proposal_index = 0;
    size_t simulate_unadjusted_proposal_index = 0;
    size_t simulate_imh_proposal_index = 0;
    size_t simulate_m_proposal_index = 0;

    for (auto i = order_of_mcmc.begin(); i != order_of_mcmc.end(); ++i) {

      if ((*i) == 1) {
        ProposalKernel *proposal;

        if (model.contains("mh_proposal")) {
          RecipeValue proposal_infos = model["mh_proposal"];
          if ((proposal_infos[simulate_mh_proposal_index]).is_list()) {
            proposal =
                get_mh_proposal(proposal_infos[simulate_mh_proposal_index],
                                model_parameters, data);
            attach_transform(proposal, proposal_infos[simulate_mh_proposal_index]);
            attach_mcmc_adaptor(proposal,
                                proposal_infos[simulate_mh_proposal_index],
                                "mh", model_parameters);
            attach_smc_adaptor(proposal, proposal_infos[simulate_mh_proposal_index],
                               model_parameters);
          } else {
            stop("Error in mcmc proposals part of model file.");
          }

          // the No-U-Turn sampler: its own move, with the Hamiltonian system
          // of its proposal (C210)
          RecipeValue mh_info = proposal_infos[simulate_mh_proposal_index];
          bool nuts = mh_info.contains("mh_proposal") && mh_info["mh_proposal"].is_list() &&
                      mh_info["mh_proposal"].contains("type") && (mh_info["mh_proposal"]["type"].as_string() == "nuts");
          if (nuts) {
            HMCProposalKernel *system = dynamic_cast<HMCProposalKernel*>(proposal);
            if (order_of_mcmc.size() > 1)
              mcmc = new NoUTurnMCMC(1, system);
            else
              mcmc = new NoUTurnMCMC(termination_method, system);
          } else if (order_of_mcmc.size() > 1) {
            mcmc = new MetropolisHastingsMCMC(1, proposal);
          } else {
            mcmc = new MetropolisHastingsMCMC(termination_method, proposal);
          }

          // The factor indices are stored inside this proposal's block of the
          // recipe (model$mh_proposal[[k]]$mh_factor_index), not in the list.
          RecipeValue current_proposal_info = proposal_infos[simulate_mh_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "mh", full_index,
              check_factor_index, simulate_mh_proposal_index + 1));

          simulate_mh_proposal_index = simulate_mh_proposal_index + 1;
        } else {
          stop("No mh_proposal found.");
        }
      } else if ((*i) == 2) {
        IndependentProposalKernel *proposal;

        if (model.contains("imh_proposal")) {
          RecipeValue proposal_infos = model["imh_proposal"];

          if ((proposal_infos[simulate_imh_proposal_index]).is_list()) {
            proposal =
                get_imh_proposal(proposal_infos[simulate_imh_proposal_index],
                                 model_parameters, data);
            attach_smc_adaptor(proposal, proposal_infos[simulate_imh_proposal_index],
                               model_parameters);
          } else {
            stop("Error in mcmc proposals part of model file.");
          }

          if (order_of_mcmc.size() > 1) {
            mcmc = new MetropolisHastingsMCMC(1, proposal);
          } else {
            mcmc = new MetropolisHastingsMCMC(termination_method, proposal);
          }

          // The factor indices are stored inside this proposal's block of the
          // recipe (model$imh_proposal[[k]]$imh_factor_index), not in the list.
          RecipeValue current_proposal_info = proposal_infos[simulate_imh_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "imh", full_index,
              check_factor_index, simulate_imh_proposal_index + 1));

          simulate_imh_proposal_index = simulate_imh_proposal_index + 1;
        } else {
          stop("No imh_proposal found.");
        }
      } else if ((*i) == 3) {
        SymmetricProposalKernel *proposal;

        if (model.contains("m_proposal")) {
          RecipeValue proposal_infos = model["m_proposal"];

          if ((proposal_infos[simulate_m_proposal_index]).is_list() &&
              proposal_infos[simulate_m_proposal_index].contains("m_elements")) {
            // a scan over elements (I4c): one Metropolis move per element, in
            // order, each with the factor indices of the move's entry
            RecipeValue current_proposal_info = proposal_infos[simulate_m_proposal_index];
            const RecipeValue &elements = current_proposal_info["m_elements"];
            std::string variable = elements["variable"].as_string();
            std::string element_type = elements["type"].as_string();
            ElementStep step = (element_type == "norm_rw") ? ElementStep::normal
                             : (element_type == "unif_rw") ? ElementStep::uniform
                             : (element_type == "gibbs_normal") ? ElementStep::normal_conditional : ElementStep::flip;
            // gibbs: over the element's states, or from its normal full
            // conditional (gibbs_normal, C206); no move through proposed particles
            bool gibbs = (element_type == "gibbs") || (element_type == "gibbs_normal");
            std::vector<double> states;
            if (element_type == "gibbs")
              states = elements["values"].as_doubles();
            std::vector<double> rows = elements["rows"].as_doubles();
            std::vector<double> columns = elements["columns"].as_doubles();
            std::vector<double> scales = elements["scales"].as_doubles();
            std::vector<MCMC *> element_moves;
            element_moves.reserve(rows.size());
            MoveFactors first_factors;
            // (gibbs has no move through proposed particles)
            // adaptation (dev/plan.md, "adaptation of element moves"): a step per
            // element adapted along the chain (adapt = scale(…)), and/or set from
            // SMC particles (smc_adapt = gaussian(c))
            bool element_adapts = elements.contains("adapt");
            double element_target = 0.44, element_kappa = 0.6, element_gain = 1.0;
            bool element_averaged = false;
            if (element_adapts) {
              BuiltinArguments adaptor = bind_builtin(elements["adapt"], BuiltinRole::MCMCAdaptor, model_parameters);
              if (adaptor.builtin() != "scale")
                stop("moves on elements adapt with adapt = scale(…).");
              if (adaptor.given("target"))
                element_target = adaptor.real("target");
              element_kappa = adaptor.real("kappa");
              element_gain = adaptor.real("C");
              element_averaged = adaptor.logical("averaged");
            }
            double element_smc_constant = 0.0;
            if (elements.contains("smc_adapt")) {
              BuiltinArguments adaptor = bind_builtin(elements["smc_adapt"], BuiltinRole::SMCAdaptor, model_parameters);
              if (adaptor.builtin() != "gaussian")
                stop("moves on elements are set from SMC particles with smc_adapt = gaussian(c).");
              element_smc_constant = adaptor.given("c") ? adaptor.real("c") : 2.38;
            }
            for (size_t e = 0; (e < rows.size()) && !gibbs; ++e) {
              ElementRandomWalkProposalKernel *element_kernel =
                  new ElementRandomWalkProposalKernel(variable, size_t(rows[e]), size_t(columns[e]), scales[e], step);
              if (element_adapts)
              {
                MCMCAdaptor* element_adaptor = new ElementScaleAdaptor(element_target, element_kappa, element_gain,
                                                                       element_averaged);
                if (elements.contains("adapt_until"))
                  element_adaptor->set_until(AdaptationCriterion::from_recipe(elements["adapt_until"]));
                element_kernel->set_mcmc_adaptor(element_adaptor);
              }
              if (element_smc_constant > 0.0)
                element_kernel->set_smc_adaptor(new ElementSMCAdaptor(element_smc_constant));
              MCMC *element_move = new MetropolisMCMC(1, element_kernel);
              MoveFactors factors = set_move_factor_indices(
                  element_move, current_proposal_info, "m", full_index,
                  check_factor_index, simulate_m_proposal_index + 1);
              if (e == 0)
                first_factors = factors;
              element_moves.push_back(element_move);
            }
            // in place (I4c-2), with the moves above for where that does not
            // apply; or (internal option, for tests) the moves above alone
            bool in_place = gibbs || !(model.contains("in_place_moves") && !model["in_place_moves"].as_bool());
            MCMCTermination *scan_termination = (order_of_mcmc.size() > 1) ? NULL : termination_method;
            if (in_place) {
              std::vector<size_t> element_rows, element_columns;
              for (size_t e = 0; e < rows.size(); ++e) {
                element_rows.push_back(size_t(rows[e]));
                element_columns.push_back(size_t(columns[e]));
              }
              InPlaceElementScanMCMC *scan = new InPlaceElementScanMCMC(scan_termination, element_moves, variable,
                                                                        element_rows, element_columns, scales, step, states);
              scan->set_smc_constant(element_smc_constant);
              mcmc = scan;
              move_factors.push_back(set_move_factor_indices(mcmc, current_proposal_info, "m", full_index,
                                                             check_factor_index, simulate_m_proposal_index + 1));
            } else {
              move_factors.push_back(first_factors);
              if (scan_termination == NULL)
                mcmc = new DeterministicScanMCMC(element_moves);
              else
                mcmc = new DeterministicScanMCMC(termination_method, element_moves);
            }
            simulate_m_proposal_index = simulate_m_proposal_index + 1;
            moves.push_back(mcmc);
            continue;
          }
          if ((proposal_infos[simulate_m_proposal_index]).is_list()) {
            proposal = get_m_proposal(proposal_infos[simulate_m_proposal_index],
                                      model_parameters, data);
            attach_transform(proposal, proposal_infos[simulate_m_proposal_index]);
            attach_mcmc_adaptor(proposal,
                                proposal_infos[simulate_m_proposal_index],
                                "m", model_parameters);
            attach_smc_adaptor(proposal, proposal_infos[simulate_m_proposal_index],
                               model_parameters);
          } else {
            stop("Error in mcmc proposals part of model file.");
          }

          if (order_of_mcmc.size() > 1) {
            mcmc = new MetropolisMCMC(1, proposal);
          } else {
            mcmc = new MetropolisMCMC(termination_method, proposal);
          }

          // The factor indices are stored inside this proposal's block of the
          // recipe (model$m_proposal[[k]]$m_factor_index), not in the list.
          RecipeValue current_proposal_info = proposal_infos[simulate_m_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "m", full_index,
              check_factor_index, simulate_m_proposal_index + 1));

          simulate_m_proposal_index = simulate_m_proposal_index + 1;
        } else {
          stop("No m_proposal found.");
        }
      } else if ((*i) == 4) {
        ProposalKernel *proposal;

        if (model.contains("unadjusted_proposal")) {
          RecipeValue proposal_infos = model["unadjusted_proposal"];
          if ((proposal_infos[simulate_unadjusted_proposal_index]).is_list()) {
            proposal = get_unadjusted_proposal(
                proposal_infos[simulate_unadjusted_proposal_index],
                model_parameters, data);
            attach_transform(proposal, proposal_infos[simulate_unadjusted_proposal_index]);
            attach_smc_adaptor(proposal,
                               proposal_infos[simulate_unadjusted_proposal_index],
                               model_parameters);
          } else {
            stop("Error in mcmc proposals part of model file.");
          }

          if (order_of_mcmc.size() > 1) {
            mcmc = new UnadjustedMCMC(1, proposal);
          } else {
            mcmc = new UnadjustedMCMC(termination_method, proposal);
          }

          // The factor indices are stored inside this proposal's block of the
          // recipe (model$unadjusted_proposal[[k]]$unadjusted_factor_index), not in the list.
          RecipeValue current_proposal_info = proposal_infos[simulate_unadjusted_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "unadjusted", full_index,
              check_factor_index, simulate_unadjusted_proposal_index + 1));

          simulate_unadjusted_proposal_index =
              simulate_unadjusted_proposal_index + 1;
        } else {
          stop("No unadjusted_proposal found.");
        }
      } else {
        Rcpp::stop("get_mcmc: invalid type in order_of_mcmc.");
      }

      moves.push_back(mcmc);
    }

    check_unadjusted_compatibility(move_factors);

    if (mcmc_weights_method.contains("func")) {
      std::vector<double> mcmc_weights;

      RecipeValue mcmc_weights_value = mcmc_weights_method["func"];
      mcmc_weights = Rcpp::as<std::vector<double>>((mcmc_weights_value).function_pointer<NumericVectorPtr>()());

      if (moves.size() == 1) {
        mcmc->set_index_if_null(full_index);
        mcmc->set_thin(mcmc_thin);
        return mcmc;
      } else {
        MCMC *new_mcmc =
            new StochasticScanMCMC(termination_method, moves, mcmc_weights);
        new_mcmc->set_index_if_null(full_index);
        new_mcmc->set_thin(mcmc_thin);
        return new_mcmc;
      }
    } else {
      if (moves.size() == 1) {
        mcmc->set_index_if_null(full_index);
        mcmc->set_thin(mcmc_thin);
        return mcmc;
      } else {
        MCMC *new_mcmc = new DeterministicScanMCMC(termination_method, moves);
        new_mcmc->set_index_if_null(full_index);
        new_mcmc->set_thin(mcmc_thin);
        return new_mcmc;
      }
    }
  } else {
    Rcpp::stop("get_mcmc: order_of_mcmc not specified in model.");
  }
}



// [[Rcpp::export]]
size_t ilike_rdtsc() { return rdtsc(); }



SMCCriterion *get_resampling_method(const RecipeValue &model,
                                    const RecipeValue &model_parameters,
                                    size_t number_of_particles) {
  if (model.contains("method")) {
    RecipeValue method_info = model["method"];
    for (size_t i = 0; i < method_info.size(); ++i) {
      RecipeValue current_method = method_info[i];
      if (current_method.contains("adaptive_resampling"))
        return make_builtin_smc_criterion(current_method["adaptive_resampling"],
                                          BuiltinRole::AdaptiveResampling,
                                          model_parameters, number_of_particles);
    }
  }

  inform_once("No method set for adaptive resampling: defaulting to resampling whenever the ESS drops below the number of particles.");
  return new ESSSMCCriterion(double(number_of_particles));
}

SMCCriterion *get_adaptive_target_method(const RecipeValue &model,
                                         const RecipeValue &model_parameters,
                                         const RecipeValue &adaptive_target_method,
                                         size_t number_of_particles) {
  if (!adaptive_target_method.contains("method"))
    stop("No valid method found for adaptive target (need method).");
  return make_builtin_smc_criterion(adaptive_target_method, BuiltinRole::AdaptiveTarget,
                                    model_parameters, number_of_particles);
}

SMCTermination *get_smc_termination_method(const RecipeValue &model,
                                           const RecipeValue &model_parameters,
                                           const RecipeValue &smc_termination_method) {
  if (smc_termination_method.length() == 0)
    return NULL;
  if (!smc_termination_method.contains("method"))
    stop("No valid method found for SMC termination (method needed).");
  return make_builtin_smc_termination(smc_termination_method, model_parameters);
}

EnsembleShifter *get_enk_shifter_method(const RecipeValue &model,
                                        const RecipeValue &model_parameters) {
  if (model.contains("method")) {
    RecipeValue method_info = model["method"];
    for (size_t i = 0; i < method_info.size(); ++i) {
      RecipeValue current_method = method_info[i];
      if (current_method.contains("enk_shifter") &&
          current_method["enk_shifter"].contains("type"))
        return make_builtin_enk_shifter(current_method["enk_shifter"], model_parameters);
    }
  }

  inform_once("No method set for shifting ensemble in EnK: defaulting to stochastic approach.");
  return new StochasticEnsembleShifter();
}

// The SMC sequences (/***smc_sequence,ilike::annealing(power,c(0,0.5,1))***/;
// the registry's role SMCSequence). The R binding evaluates the schedule and
// factor expressions (they are R expressions) and passes one built-in call
// per sequence, {builtin, arguments}, under "sequences" (role: SMCSequence,
// or EnKSequence for EnKI, which allows only annealing). Returned as before:
// (types, variables, schedules, bisections[, factors affected]).
RecipeValue get_smc_sequencer_info(const RecipeValue &model, const RecipeValue &model_parameters,
                            const RecipeValue &smc_sequencer_method, BuiltinRole role) {
  std::vector<std::string> types;
  std::vector<std::string> variables;
  std::vector<std::vector<double>> schedules;
  std::vector<double> factors_affected;
  bool any_factors = false;
  size_t number_of_bisections = 100;

  if (smc_sequencer_method.contains("sequences")) {
    const RecipeValue &sequences = smc_sequencer_method["sequences"];
    for (size_t i = 0; i < sequences.size(); ++i) {
      BuiltinArguments arguments =
          bind_builtin(sequences[i], role, model_parameters);
      types.push_back(arguments.builtin());
      variables.push_back(arguments.variable("variable"));
      arma::colvec schedule = arguments.vector("schedule");
      schedules.push_back(std::vector<double>(schedule.begin(), schedule.end()));
      if (arguments.given("factors")) {
        any_factors = true;
        arma::colvec factors = arguments.vector("factors");
        factors_affected.insert(factors_affected.end(), factors.begin(), factors.end());
      }
    }
    if (!sequences.is_null() && sequences.size() > 0)
      inform_once("Bisections not set in adaptive target method, using 100.");
  }

  if (any_factors)
    return RecipeValue::list({types, variables, schedules, number_of_bisections,
                              factors_affected});
  return RecipeValue::list({types, variables, schedules, number_of_bisections});
}
