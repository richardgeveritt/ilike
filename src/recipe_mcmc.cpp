#include "recipe_internal.h"

IndependentProposalKernel *
get_proposal(const List &model, const List &model_parameters, Data *data) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.containsElementNamed("is_proposal")) {
    List is_proposals = model["is_proposal"];

    for (size_t i = 0; i < is_proposals.size(); ++i) {
      if (Rf_isNewList(is_proposals[i])) {
        List current_proposal = is_proposals[i];
        if (current_proposal.containsElementNamed("is_proposal")) {
          List proposal_info = current_proposal["is_proposal"];

          if (Rf_isNewList(proposal_info)) {
            if (proposal_info.containsElementNamed("type") &&
                proposal_info.containsElementNamed("variables")) {
              std::string type = proposal_info["type"];
              if (iequals(type, "norm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, proposal_info, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, sd);
              } else if (iequals(type, "mvnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, proposal_info, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, cov);
              } else if (iequals(type, "lnorm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, proposal_info, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, sd);
              } else if (iequals(type, "mvlnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, proposal_info, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, cov);
              } else if (iequals(type, "gamma")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, proposal_info, type);

                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];

                proposal =
                    new GammaIndependentProposalKernel(variable, shape, rate);
              } else {
                Rcout << "Proposal type " << type;
                stop("Proposal type unknown");
              }
            } else {
              stop("Missing information for proposal in model file.");
            }
          } else {
            stop("Error in proposal section of model file.");
          }
        } else if (current_proposal.containsElementNamed(
                       "evaluate_log_is_proposal") &&
                   current_proposal.containsElementNamed(
                       "simulate_is_proposal") &&
                   current_proposal.containsElementNamed("type")) {
          SEXP evaluate_log_is_proposal_SEXP =
              current_proposal["evaluate_log_is_proposal"];
          SEXP simulate_is_proposal_SEXP =
              current_proposal["simulate_is_proposal"];
          size_t type = current_proposal["type"];
          if (type == 1) {
            proposal = new CustomDistributionProposalKernel(
                load_simulate_distribution(simulate_is_proposal_SEXP),
                load_evaluate_log_distribution(evaluate_log_is_proposal_SEXP));
          } else if (type == 2) {
            proposal = new CustomIndependentProposalKernel(
                load_simulate_independent_proposal(simulate_is_proposal_SEXP),
                load_evaluate_log_independent_proposal(
                    evaluate_log_is_proposal_SEXP));
          } else if (type == 3) {
            proposal = new CustomGuidedDistributionProposalKernel(
                load_simulate_guided_distribution(simulate_is_proposal_SEXP),
                load_evaluate_log_guided_distribution(
                    evaluate_log_is_proposal_SEXP),
                data);
          } else if (type == 4) {
            proposal = new CustomGuidedIndependentProposalKernel(
                load_simulate_guided_independent_proposal(
                    simulate_is_proposal_SEXP),
                load_evaluate_log_guided_independent_proposal(
                    evaluate_log_is_proposal_SEXP),
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

IndependentProposalKernel *get_prior_as_proposal(const List &model,
                                                 const List &model_parameters) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.containsElementNamed("factor")) {
    List factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {
      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];
        if (current_factor.containsElementNamed("prior")) {
          if (Rf_isNewList(current_factor["prior"])) {
            List current_prior = current_factor["prior"];
            if (current_prior.containsElementNamed("type") &&
                current_prior.containsElementNamed("variables")) {
              std::string type = current_prior["type"];

              if (iequals(type, "norm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, sd);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "mvnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, cov);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "lnorm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, sd);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "mvlnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, cov);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "gamma")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];

                proposal =
                    new GammaIndependentProposalKernel(variable, shape, rate);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else {
                Rcout << "Prior type " << type;
                stop("Prior type unknown");
              }
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if ((current_factor.containsElementNamed(
                       "evaluate_log_prior")) &&
                   (current_factor.containsElementNamed("simulate_prior"))) {
          SEXP evaluate_log_prior_SEXP = current_factor["evaluate_log_prior"];

          SEXP simulate_prior_SEXP = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(
              load_simulate_distribution(simulate_prior_SEXP),
              load_evaluate_log_distribution(evaluate_log_prior_SEXP));

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
get_prior_as_simulate_only_proposal(const List &model,
                                    const List &model_parameters) {
  IndependentProposalKernel *proposal = NULL;
  std::vector<IndependentProposalKernel *> proposals;

  if (model.containsElementNamed("factor")) {
    List factors = model["factor"];

    for (size_t i = 0; i < factors.size(); ++i) {
      if (Rf_isNewList(factors[i])) {
        List current_factor = factors[i];
        if (current_factor.containsElementNamed("prior")) {
          if (Rf_isNewList(current_factor["prior"])) {
            List current_prior = current_factor["prior"];
            if (current_prior.containsElementNamed("type") &&
                current_prior.containsElementNamed("variables")) {
              std::string type = current_prior["type"];

              if (iequals(type, "norm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, sd);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "mvnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal =
                    new GaussianIndependentProposalKernel(variable, mean, cov);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "lnorm")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double mean = info[1];
                double sd = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, sd);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "mvlnorm")) {
                List info =
                    get_single_variable_vector_and_matrix_parameter_info(
                        model_parameters, current_prior, type);

                std::string variable = info[0];
                arma::colvec mean = info[1];
                arma::mat cov = info[2];

                proposal = new LogGaussianIndependentProposalKernel(variable,
                                                                    mean, cov);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else if (iequals(type, "gamma")) {
                List info = get_single_variable_two_double_parameter_info(
                    model_parameters, current_prior, type);

                std::string variable = info[0];
                double shape = info[1];
                double rate = info[2];

                proposal =
                    new GammaIndependentProposalKernel(variable, shape, rate);

                if (proposal != NULL) {
                  proposals.push_back(proposal);
                }
              } else {
                Rcout << "Prior type " << type;
                stop("Prior type unknown");
              }
            } else {
              stop("Missing information for prior in model file.");
            }
          } else {
            stop("Error in factors section of model file.");
          }
        } else if (current_factor.containsElementNamed("simulate_prior")) {
          SEXP simulate_prior_SEXP = current_factor["simulate_prior"];
          proposal = new CustomDistributionProposalKernel(
              load_simulate_distribution(simulate_prior_SEXP));

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

ProposalKernel *get_mh_proposal(const List &current_proposal,
                                const List &model_parameters, Data *data) {
  ProposalKernel *proposal;
  if (current_proposal.containsElementNamed("mh_proposal")) {
    List proposal_info = current_proposal["mh_proposal"];

    if (Rf_isNewList(proposal_info)) {
      if (proposal_info.containsElementNamed("type") &&
          proposal_info.containsElementNamed("variables")) {
        std::string type = proposal_info["type"];
        if (iequals(type, "norm_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new GaussianRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "mvnorm_rw")) {
          // (variable, Sigma) or (variable, Sigma, s^2); s^2 defaults to 1.
          List info =
              get_single_variable_matrix_and_optional_double_parameter_info(
                  model_parameters, proposal_info, type, 1.0);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new GaussianRandomWalkProposalKernel(variable, cov, scale);
        } else if (iequals(type, "unif_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "mvunif_rw")) {
          List info = get_single_variable_vector_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec halfwidth = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, halfwidth);
        } else if (iequals(type, "langevin")) {
          List info = get_single_variable_matrix_and_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new LangevinProposalKernel(variable, cov, scale,
                                                new DirectGradientEstimator());
        } else if (iequals(type, "hamiltonian")) {
          // HMCProposalKernel is a stub (simulate returns no parameters, so
          // the chain never moves); refuse it rather than return a stuck chain.
          stop("MH proposal ilike::hamiltonian (HMC) is not written yet.");
          List info = get_single_variable_matrix_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];

          proposal = new HMCProposalKernel(variable, cov,
                                           new DirectGradientEstimator());
        } else if (iequals(type, "barker")) {
          List info = get_single_variable_matrix_and_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new BarkerDynamicsProposalKernel(
              variable, cov, scale, new DirectGradientEstimator());
        } else if (iequals(type, "mirror")) {
          List info =
              get_single_variable_vector_and_matrix_and_double_parameter_info(
                  model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          double scale = info[3];

          proposal = new MirrorProposalKernel(variable, mean, cov, scale);
        } else {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.containsElementNamed(
                 "evaluate_log_mh_proposal") &&
             current_proposal.containsElementNamed("simulate_mh_proposal") &&
             current_proposal.containsElementNamed("type")) {
    SEXP evaluate_log_mh_proposal_SEXP =
        current_proposal["evaluate_log_mh_proposal"];
    SEXP simulate_mh_proposal_SEXP = current_proposal["simulate_mh_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsProposalKernel(
          load_simulate_no_params_mcmc_proposal(simulate_mh_proposal_SEXP),
          load_evaluate_log_no_params_mcmc_proposal(
              evaluate_log_mh_proposal_SEXP));
    } else if (type == 2) {
      proposal = new CustomProposalKernel(
          load_simulate_mcmc_proposal(simulate_mh_proposal_SEXP),
          load_evaluate_log_mcmc_proposal(evaluate_log_mh_proposal_SEXP));
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsProposalKernel(
          load_simulate_guided_no_params_mcmc_proposal(
              simulate_mh_proposal_SEXP),
          load_evaluate_log_guided_no_params_mcmc_proposal(
              evaluate_log_mh_proposal_SEXP),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedProposalKernel(
          load_simulate_guided_mcmc_proposal(simulate_mh_proposal_SEXP),
          load_evaluate_log_guided_mcmc_proposal(evaluate_log_mh_proposal_SEXP),
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

ProposalKernel *get_unadjusted_proposal(const List &current_proposal,
                                        const List &model_parameters,
                                        Data *data) {
  ProposalKernel *proposal;
  if (current_proposal.containsElementNamed("unadjusted_proposal")) {
    List proposal_info = current_proposal["unadjusted_proposal"];

    if (Rf_isNewList(proposal_info)) {
      if (proposal_info.containsElementNamed("type") &&
          proposal_info.containsElementNamed("variables")) {
        std::string type = proposal_info["type"];
        if (iequals(type, "norm_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new GaussianRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "mvnorm_rw")) {
          // (variable, Sigma) or (variable, Sigma, s^2); s^2 defaults to 1.
          List info =
              get_single_variable_matrix_and_optional_double_parameter_info(
                  model_parameters, proposal_info, type, 1.0);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new GaussianRandomWalkProposalKernel(variable, cov, scale);
        } else if (iequals(type, "unif_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "mvunif_rw")) {
          List info = get_single_variable_vector_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec halfwidth = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, halfwidth);
        } else if (iequals(type, "langevin")) {
          List info = get_single_variable_matrix_and_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new LangevinProposalKernel(variable, cov, scale,
                                                new DirectGradientEstimator());
        } else if (iequals(type, "hamiltonian")) {
          // HMCProposalKernel is a stub (simulate returns no parameters, so
          // the chain never moves); refuse it rather than return a stuck chain.
          stop("MH proposal ilike::hamiltonian (HMC) is not written yet.");
          List info = get_single_variable_matrix_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];

          proposal = new HMCProposalKernel(variable, cov,
                                           new DirectGradientEstimator());
        } else if (iequals(type, "barker")) {
          List info = get_single_variable_matrix_and_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new BarkerDynamicsProposalKernel(
              variable, cov, scale, new DirectGradientEstimator());
        } else if (iequals(type, "mirror")) {
          List info =
              get_single_variable_vector_and_matrix_and_double_parameter_info(
                  model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          double scale = info[3];

          proposal = new MirrorProposalKernel(variable, mean, cov, scale);
        } else {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.containsElementNamed(
                 "simulate_unadjusted_proposal") &&
             current_proposal.containsElementNamed("type")) {
    SEXP simulate_unadjusted_proposal_SEXP =
        current_proposal["simulate_unadjusted_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsProposalKernel(
          load_simulate_no_params_mcmc_proposal(
              simulate_unadjusted_proposal_SEXP));
    } else if (type == 2) {
      proposal = new CustomProposalKernel(
          load_simulate_mcmc_proposal(simulate_unadjusted_proposal_SEXP));
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsProposalKernel(
          load_simulate_guided_no_params_mcmc_proposal(
              simulate_unadjusted_proposal_SEXP),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedProposalKernel(
          load_simulate_guided_mcmc_proposal(simulate_unadjusted_proposal_SEXP),
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

IndependentProposalKernel *get_imh_proposal(const List &current_proposal,
                                            const List &model_parameters,
                                            Data *data) {
  IndependentProposalKernel *proposal;
  if (current_proposal.containsElementNamed("imh_proposal")) {
    List proposal_info = current_proposal["imh_proposal"];

    if (Rf_isNewList(proposal_info)) {
      if (proposal_info.containsElementNamed("type") &&
          proposal_info.containsElementNamed("variables")) {
        std::string type = proposal_info["type"];
        if (iequals(type, "norm")) {
          List info = get_single_variable_two_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double mean = info[1];
          double sd = info[2];

          proposal = new GaussianIndependentProposalKernel(variable, mean, sd);
        } else if (iequals(type, "mvnorm")) {
          List info = get_single_variable_vector_and_matrix_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];

          proposal = new GaussianIndependentProposalKernel(variable, mean, cov);
        } else if (iequals(type, "lnorm")) {
          List info = get_single_variable_two_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double mean = info[1];
          double sd = info[2];

          proposal =
              new LogGaussianIndependentProposalKernel(variable, mean, sd);
        } else if (iequals(type, "mvlnorm")) {
          List info = get_single_variable_vector_and_matrix_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];

          proposal =
              new LogGaussianIndependentProposalKernel(variable, mean, cov);
        } else if (iequals(type, "gamma")) {
          List info = get_single_variable_two_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double shape = info[1];
          double rate = info[2];

          proposal = new GammaIndependentProposalKernel(variable, shape, rate);
        } else {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.containsElementNamed(
                 "evaluate_log_imh_proposal") &&
             current_proposal.containsElementNamed("simulate_imh_proposal") &&
             current_proposal.containsElementNamed("type")) {
    SEXP evaluate_log_imh_proposal_SEXP =
        current_proposal["evaluate_log_imh_proposal"];
    SEXP simulate_imh_proposal_SEXP = current_proposal["simulate_imh_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomDistributionProposalKernel(
          load_simulate_distribution(simulate_imh_proposal_SEXP),
          load_evaluate_log_distribution(evaluate_log_imh_proposal_SEXP));
    } else if (type == 2) {
      proposal = new CustomIndependentProposalKernel(
          load_simulate_independent_proposal(simulate_imh_proposal_SEXP),
          load_evaluate_log_independent_proposal(
              evaluate_log_imh_proposal_SEXP));
    } else if (type == 3) {
      proposal = new CustomGuidedDistributionProposalKernel(
          load_simulate_guided_distribution(simulate_imh_proposal_SEXP),
          load_evaluate_log_guided_distribution(evaluate_log_imh_proposal_SEXP),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedIndependentProposalKernel(
          load_simulate_guided_independent_proposal(simulate_imh_proposal_SEXP),
          load_evaluate_log_guided_independent_proposal(
              evaluate_log_imh_proposal_SEXP),
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

SymmetricProposalKernel *get_m_proposal(const List &current_proposal,
                                        const List &model_parameters,
                                        Data *data) {
  SymmetricProposalKernel *proposal;
  if (current_proposal.containsElementNamed("m_proposal")) {
    List proposal_info = current_proposal["m_proposal"];

    if (Rf_isNewList(proposal_info)) {
      if (proposal_info.containsElementNamed("type") &&
          proposal_info.containsElementNamed("variables")) {
        std::string type = proposal_info["type"];
        if (iequals(type, "norm_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new GaussianRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "pc_rw")) {
          // (variable, Sigma[, s^2[, m]]): moves along one of the first m
          // principal directions of Sigma at a time (Andrieu & Thoms 2008,
          // Algorithm 8); m = 0 or absent: all of them
          size_t number_of_directions = 0;
          List pc_info = clone(proposal_info);
          if (pc_info.containsElementNamed("parameters")) {
            List pc_parameters = pc_info["parameters"];
            if (pc_parameters.size() > 3)
              stop("ilike::pc_rw: too many arguments.");
            if (pc_parameters.size() == 3) {
              double m = extract_double_parameter(pc_parameters, model_parameters, 2);
              if ((m < 0.0) || (m != floor(m)))
                stop("ilike::pc_rw: the number of directions must be a whole number.");
              number_of_directions = size_t(m);
              pc_info["parameters"] = List::create(pc_parameters[0], pc_parameters[1]);
            }
          }
          List info =
              get_single_variable_matrix_and_optional_double_parameter_info(
                  model_parameters, pc_info, type, 1.0);
          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];
          proposal = new PrincipalComponentsProposalKernel(variable, cov, scale,
                                                           number_of_directions);
        } else if (iequals(type, "componentwise_rw")) {
          // (variable, Sigma) or (variable, Sigma, s^2): one component k at
          // a time, N(0, s^2 [Sigma]_kk) (Andrieu & Thoms 2008, Algorithm 5)
          List info =
              get_single_variable_matrix_and_optional_double_parameter_info(
                  model_parameters, proposal_info, type, 1.0);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new ComponentwiseRandomWalkProposalKernel(variable, cov, scale);
        } else if (iequals(type, "mvnorm_rw")) {
          // (variable, Sigma) or (variable, Sigma, s^2); s^2 defaults to 1.
          List info =
              get_single_variable_matrix_and_optional_double_parameter_info(
                  model_parameters, proposal_info, type, 1.0);

          std::string variable = info[0];
          arma::mat cov = info[1];
          double scale = info[2];

          proposal = new GaussianRandomWalkProposalKernel(variable, cov, scale);
        } else if (iequals(type, "unif_rw")) {
          List info = get_single_variable_one_double_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          double sd = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, sd);
        } else if (iequals(type, "mvunif_rw")) {
          List info = get_single_variable_vector_parameter_info(
              model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec halfwidth = info[1];

          proposal = new UniformRandomWalkProposalKernel(variable, halfwidth);
        } else if (iequals(type, "mirror")) {
          List info =
              get_single_variable_vector_and_matrix_and_double_parameter_info(
                  model_parameters, proposal_info, type);

          std::string variable = info[0];
          arma::colvec mean = info[1];
          arma::mat cov = info[2];
          double scale = info[3];

          proposal = new MirrorProposalKernel(variable, mean, cov, scale);
        } else {
          Rcout << "Proposal type " << type;
          stop("Proposal type unknown");
        }
      } else {
        stop("Missing information for proposal in model file.");
      }
    } else {
      stop("Error in proposal section of model file.");
    }
  } else if (current_proposal.containsElementNamed("simulate_m_proposal") &&
             current_proposal.containsElementNamed("type")) {
    SEXP simulate_m_proposal_SEXP = current_proposal["simulate_m_proposal"];
    size_t type = current_proposal["type"];
    if (type == 1) {
      proposal = new CustomNoParamsSymmetricProposalKernel(
          load_simulate_no_params_mcmc_proposal(simulate_m_proposal_SEXP));
    } else if (type == 2) {
      proposal = new CustomSymmetricProposalKernel(
          load_simulate_mcmc_proposal(simulate_m_proposal_SEXP));
    } else if (type == 3) {
      proposal = new CustomGuidedNoParamsSymmetricProposalKernel(
          load_simulate_guided_no_params_mcmc_proposal(
              simulate_m_proposal_SEXP),
          data);
    } else if (type == 4) {
      proposal = new CustomGuidedSymmetricProposalKernel(
          load_simulate_guided_mcmc_proposal(simulate_m_proposal_SEXP), data);
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

MCMCTermination *make_mcmc_termination(const List &model_parameters,
                                       const List &mcmc_termination_method) {
  MCMCTermination *termination_method;

  if (mcmc_termination_method.containsElementNamed("method") &&
      mcmc_termination_method.containsElementNamed("values")) {
    std::string method = mcmc_termination_method["method"];
    if (iequals(method, "iterations")) {
      if (Rf_isNewList(mcmc_termination_method["values"])) {
        List values = mcmc_termination_method["values"];

        if (values.length() == 1) {
          int number_of_iterations =
              extract_int_parameter(values, model_parameters, 0);
          termination_method =
              new IterationsMCMCTermination(number_of_iterations);
        } else {
          stop("MCMC termination using iterations requires you to specify only "
               "the number of iterations.");
        }
      } else {
        stop("MCMC termination using iterations requires you to specify only "
             "the number of iterations.");
      }
    } else if (iequals(method, "se")) {
      // SEMCMCTermination::terminate() ignores its threshold and stops only
      // at the maximum number of iterations (never, if no maximum is given),
      // so refuse it until the standard-error rule is written. The parsing
      // below is kept for that.
      stop("MCMC termination ilike::se is not written yet: use "
           "ilike::iterations.");
      if (Rf_isNewList(mcmc_termination_method["values"])) {
        List values = mcmc_termination_method["values"];

        if (values.length() == 1) {
          double threshold =
              extract_double_parameter(values, model_parameters, 0);
          int number_of_iterations = std::numeric_limits<int>::max();

          termination_method =
              new SEMCMCTermination(threshold, number_of_iterations);
        } else if (values.length() == 2) {
          double threshold =
              extract_double_parameter(values, model_parameters, 0);
          int number_of_iterations =
              extract_int_parameter(values, model_parameters, 1);

          termination_method =
              new SEMCMCTermination(threshold, number_of_iterations);
        } else {
          stop("MCMC termination using standard error requires you to specify "
               "a threshold.");
        }
      } else {
        stop("MCMC termination using standard error requires you to specify a "
             "threshold.");
      }
    } else {
      stop("Invalid method for MCMC termination.");
    }
  } else {
    stop("No valid method found for MCMC termination.");
  }

  return termination_method;
}

static MoveFactors set_move_factor_indices(MCMC *mcmc, const List &proposal_info,
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
  bool has_used = proposal_info.containsElementNamed(used_name.c_str());
  bool has_cut = proposal_info.containsElementNamed(cut_name.c_str());
  if (!has_used && !has_cut)
    return move_factors;
  mcmc->set_check_factor_index(check_factor_index, move_factors.label);
  if (full_index == NULL)
    Rcpp::stop("make_mcmc: full index of factors not available.");

  std::vector<size_t> all_factors(full_index->begin(), full_index->end());

  auto read_index = [&](const std::string &name) {
    SEXP index_SEXP = proposal_info[name];
    NumericVector index = load_numeric_vector(index_SEXP);
    std::vector<size_t> result;
    for (R_xlen_t k = 0; k < index.length(); ++k) {
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
// Either a built-in adaptor block after the proposal
// (, stored in the proposal's
// part of the recipe under "<prefix>_adaptor"), or the user's adaptation
// function (, under "adapt_<prefix>_proposal").
// Built-in adaptors; the arguments after the first are optional:
//   scaled_am(theta, target, kappa, C, r)        Algorithm 4
//   am(theta, kappa, C, r)                       Algorithm 2
//   rb_am(theta, kappa, C, r)                    Algorithm 3
//   scale(theta, target, kappa, C)               the scale recursion of Algorithm 4
//   componentwise_am(theta, target, kappa, C, r) Algorithm 5 (proposal componentwise_rw)
//   pc_am(theta, target, kappa, C, r)            Algorithm 8, eigen-decomposition
//                                                every r iterations (proposal pc_rw)
//   online_pc_am(theta, target, kappa, C, r)     Algorithm 8, the online PCA
//                                                recursion every r iterations
//   acceptance_rate(name, target, kappa, C)      a named proposal parameter of
//                                                a proposal written by the user
// target: acceptance rate (default 0.234 for random walks and named
// parameters, 0.44 for componentwise and principal-components moves, for the
// mirror proposal a rate depending on the dimension (0.50 for d = 1 down to
// 0.23 for d >= 50, found by experiment), 0.574 for Langevin, 0.4 for Barker); gain
// C (n+1)^(-kappa) (defaults kappa = 0.6, C = 1); r: refresh interval
// (default 0: every d iterations; every iteration for online_pc_am). The
// mirror proposal (ilike::mirror) takes scaled_am, am, rb_am and scale, and
// its mean is set to the running estimate of the target's mean.
// Averaging (Polyak & Juditsky 1992): the argument "averaged", anywhere
// after the first, makes the kernel use the averages of the recursion's
// iterates; for a user's adaptation function, the block
//does the same.
void attach_mcmc_adaptor(ProposalKernel *proposal, const List &proposal_info,
                         const std::string &prefix,
                         const List &model_parameters) {
  std::string adaptor_name = prefix + "_adaptor";
  std::string function_name = "adapt_" + prefix + "_proposal";
  bool has_block = proposal_info.containsElementNamed(adaptor_name.c_str());
  bool has_function = proposal_info.containsElementNamed(function_name.c_str());
  if (has_function) {
    // the only block allowed with a user's function: ilike::averaged()
    bool function_averaging = false;
    if (has_block) {
      List block_info = proposal_info[adaptor_name];
      std::string block_type = block_info["type"];
      if (!iequals(block_type, "averaged"))
        stop("A proposal with an " + function_name + " function can have only "
             "the block " + adaptor_name + ",ilike::averaged() with it.");
      function_averaging = true;
    }
    SEXP adapt_SEXP = proposal_info[function_name];
    proposal->set_mcmc_adaptor(new ProposalParameterAdaptor(
        load_adapt_proposal(adapt_SEXP), 0.6, 1.0, function_averaging));
    return;
  }
  if (!has_block)
    return;

  std::string proposal_type = "";
  std::string proposal_name = prefix + "_proposal";
  if (proposal_info.containsElementNamed(proposal_name.c_str())) {
    SEXP proposal_SEXP = proposal_info[proposal_name];
    if (Rf_isNewList(proposal_SEXP)) {
      List proposal_list(proposal_SEXP);
      if (proposal_list.containsElementNamed("type"))
        proposal_type = as<std::string>(proposal_list["type"]);
    }
  }

  List adaptor_info = proposal_info[adaptor_name];
  std::string type = adaptor_info["type"];
  if (iequals(type, "averaged"))
    stop(adaptor_name + ",ilike::averaged() goes with an " + function_name +
         " function; for built-in adaptors, give \"averaged\" as an argument "
         "(e.g. ilike::scaled_am(theta,averaged)).");
  std::string variable = adaptor_info["variables"];
  List given_parameters;
  if (adaptor_info.containsElementNamed("parameters"))
    given_parameters = adaptor_info["parameters"];
  // the keyword "averaged" (Polyak-Juditsky averaging) may be anywhere
  bool averaging = false;
  List parameters;
  for (R_xlen_t a = 0; a < given_parameters.size(); ++a) {
    std::string argument = as<std::string>(given_parameters[a]);
    if (iequals(argument, "averaged"))
      averaging = true;
    else
      parameters.push_back(given_parameters[a]);
  }

  bool random_walk = iequals(proposal_type, "norm_rw") ||
                     iequals(proposal_type, "mvnorm_rw");
  bool componentwise = iequals(proposal_type, "componentwise_rw");
  bool principal_components = iequals(proposal_type, "pc_rw");
  bool mirror = iequals(proposal_type, "mirror");
  bool gradient = iequals(proposal_type, "langevin") ||
                  iequals(proposal_type, "barker");
  bool user_written = (proposal_type == "");

  double default_target = 0.234;
  if (componentwise || principal_components)
    default_target = 0.44;
  else if (mirror)
    default_target = 0.0; // by the variable's dimension (mirror_default_target_acceptance)
  else if (iequals(proposal_type, "langevin"))
    default_target = 0.574;
  else if (iequals(proposal_type, "barker"))
    default_target = 0.4;

  // the methods each proposal takes
  std::string fits;
  bool fit = false;
  if (iequals(type, "scaled_am") || iequals(type, "am") || iequals(type, "rb_am")) {
    fit = random_walk || mirror;
    fits = "the proposals ilike::norm_rw, ilike::mvnorm_rw and ilike::mirror";
  } else if (iequals(type, "scale")) {
    fit = random_walk || mirror || gradient;
    fits = "the proposals ilike::norm_rw, ilike::mvnorm_rw, ilike::mirror, ilike::langevin and ilike::barker";
  } else if (iequals(type, "componentwise_am")) {
    fit = componentwise;
    fits = "the proposal ilike::componentwise_rw";
  } else if (iequals(type, "pc_am") || iequals(type, "online_pc_am")) {
    fit = principal_components;
    fits = "the proposal ilike::pc_rw";
  } else if (iequals(type, "acceptance_rate")) {
    fit = user_written;
    fits = "proposals written by the user (with proposal_parameters)";
  } else {
    stop("Unknown MCMC adaptor ilike::" + type +
         " (available: ilike::scaled_am, ilike::am, ilike::rb_am, ilike::scale, "
         "ilike::componentwise_am, ilike::pc_am, ilike::online_pc_am, "
         "ilike::acceptance_rate).");
  }
  if (!fit) {
    if (componentwise)
      stop(adaptor_name + ",ilike::" + type + ": the proposal "
           "ilike::componentwise_rw is adapted by ilike::componentwise_am.");
    if (principal_components)
      stop(adaptor_name + ",ilike::" + type + ": the proposal "
           "ilike::pc_rw is adapted by ilike::pc_am or ilike::online_pc_am.");
    stop(adaptor_name + ",ilike::" + type + ": applies to " + fits + ".");
  }

  bool has_target = !(iequals(type, "am") || iequals(type, "rb_am"));
  bool has_refresh = !(iequals(type, "scale") || iequals(type, "acceptance_rate"));
  size_t maximum_arguments = (has_target ? 1 : 0) + 2 + (has_refresh ? 1 : 0);
  if (size_t(parameters.size()) > maximum_arguments)
    stop(adaptor_name + ",ilike::" + type + ": too many arguments.");

  size_t k = 0;
  double target = default_target;
  if (has_target) {
    if (size_t(parameters.size()) > k) {
      target = extract_double_parameter(parameters, model_parameters, k);
      if ((target <= 0.0) || (target >= 1.0))
        stop(adaptor_name + ": the target acceptance rate must be between 0 and 1.");
    }
    k = k + 1;
  }
  double kappa = 0.6;
  double constant = 1.0;
  double refresh = 0.0;
  if (size_t(parameters.size()) > k)
    kappa = extract_double_parameter(parameters, model_parameters, k);
  if (size_t(parameters.size()) > k + 1)
    constant = extract_double_parameter(parameters, model_parameters, k + 1);
  if (has_refresh && (size_t(parameters.size()) > k + 2))
    refresh = extract_double_parameter(parameters, model_parameters, k + 2);

  if (constant * pow(2.0, -kappa) >= 1.0)
    stop(adaptor_name + ": the gain C (n+1)^(-kappa) must be below 1 from the "
         "first update (C 2^(-kappa) < 1).");
  if ((refresh < 0.0) || (refresh != floor(refresh)))
    stop(adaptor_name + ": the refresh interval must be a whole number (0: the "
         "default).");

  std::vector<std::string> variables;
  variables.push_back(variable);

  if (iequals(type, "acceptance_rate")) {
    proposal->set_mcmc_adaptor(new ProposalParameterAdaptor(
        variables, target, kappa, constant, averaging));
    return;
  }

  GaussianAdaptationMethod method = GaussianAdaptationMethod::SCALED_AM;
  if (iequals(type, "am"))
    method = GaussianAdaptationMethod::AM;
  else if (iequals(type, "rb_am"))
    method = GaussianAdaptationMethod::RB_AM;
  else if (iequals(type, "scale"))
    method = GaussianAdaptationMethod::SCALE;
  else if (iequals(type, "componentwise_am"))
    method = GaussianAdaptationMethod::COMPONENTWISE_AM;
  else if (iequals(type, "pc_am"))
    method = GaussianAdaptationMethod::PC_EIGEN;
  else if (iequals(type, "online_pc_am"))
    method = GaussianAdaptationMethod::PC_ONLINE;

  proposal->set_mcmc_adaptor(new GaussianMCMCAdaptor(
      variables, method, target, kappa, constant, size_t(refresh), averaging));
}

// Adaptation of a proposal from the SMC particles (dev/plan.md 3.2, B2): the
// blockafter the proposal,
// stored in the proposal's part of the recipe under "smc_adaptor". Before
// each move of an SMC algorithm (or ensemble Kalman method), the proposal's
// covariance is set to the weighted covariance of the particles (and a
// mirror's centre to their weighted mean; a principal-components proposal's
// directions to its eigenvectors); with c, the scale to c^2/d (c^2 for
// componentwise and principal-components proposals).
void attach_smc_adaptor(ProposalKernel *proposal, const List &proposal_info,
                        const List &model_parameters) {
  if (!proposal_info.containsElementNamed("smc_adaptor"))
    return;
  List adaptor_info = proposal_info["smc_adaptor"];
  std::string type = adaptor_info["type"];
  if (!iequals(type, "gaussian"))
    stop("Unknown SMC adaptor ilike::" + type + " (available: ilike::gaussian).");
  std::string variable = adaptor_info["variables"];
  List parameters;
  if (adaptor_info.containsElementNamed("parameters"))
    parameters = adaptor_info["parameters"];
  if (parameters.size() > 1)
    stop("smc_adaptor,ilike::gaussian: too many arguments (the variable, then optionally the constant c).");
  double constant = -1.0;
  if (parameters.size() == 1) {
    constant = extract_double_parameter(parameters, model_parameters, 0);
    if (constant <= 0.0)
      stop("smc_adaptor,ilike::gaussian: the constant c must be positive.");
  }
  if (proposal->get_gaussian_proposal_info(variable) == NULL)
    stop("smc_adaptor,ilike::gaussian: the proposal has no Gaussian covariance "
         "for " + variable + " (available for ilike::norm_rw, ilike::mvnorm_rw, "
         "ilike::componentwise_rw, ilike::pc_rw, ilike::mirror, ilike::langevin, "
         "ilike::barker, and the independent proposals ilike::norm and ilike::mvnorm).");
  proposal->set_smc_adaptor(new GaussianSMCAdaptor(variable, constant));
}

MCMC *make_mcmc(const List &model, const List &model_parameters, Data *data,
                const List &mcmc_termination_method,
                const List &mcmc_weights_method, Index *full_index) {
  MCMCTermination *termination_method =
      make_mcmc_termination(model_parameters, mcmc_termination_method);

  MCMC *mcmc = NULL;

  // Optional consistency check of the factor indices, set by the R wrappers
  // (MCMC(..., check_factor_index = TRUE)) as an element of the recipe.
  bool check_factor_index = false;
  if (model.containsElementNamed("check_factor_index"))
    check_factor_index = as<bool>(model["check_factor_index"]);

  // MCMC(..., thin): store every thin-th iteration (recipe element mcmc_thin).
  size_t mcmc_thin = 1;
  if (model.containsElementNamed("mcmc_thin"))
    mcmc_thin = as<size_t>(model["mcmc_thin"]);

  if (model.containsElementNamed("order_of_mcmc")) {
    NumericVector order_of_mcmc = model["order_of_mcmc"];
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

        if (model.containsElementNamed("mh_proposal")) {
          List proposal_infos = model["mh_proposal"];
          if (Rf_isNewList(proposal_infos[simulate_mh_proposal_index])) {
            proposal =
                get_mh_proposal(proposal_infos[simulate_mh_proposal_index],
                                model_parameters, data);
            attach_mcmc_adaptor(proposal,
                                proposal_infos[simulate_mh_proposal_index],
                                "mh", model_parameters);
            attach_smc_adaptor(proposal, proposal_infos[simulate_mh_proposal_index],
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
          // recipe (model$mh_proposal[[k]]$mh_factor_index), not in the list.
          List current_proposal_info = proposal_infos[simulate_mh_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "mh", full_index,
              check_factor_index, simulate_mh_proposal_index + 1));

          simulate_mh_proposal_index = simulate_mh_proposal_index + 1;
        } else {
          stop("No mh_proposal found.");
        }
      } else if ((*i) == 2) {
        IndependentProposalKernel *proposal;

        if (model.containsElementNamed("imh_proposal")) {
          List proposal_infos = model["imh_proposal"];

          if (Rf_isNewList(proposal_infos[simulate_imh_proposal_index])) {
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
          List current_proposal_info = proposal_infos[simulate_imh_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "imh", full_index,
              check_factor_index, simulate_imh_proposal_index + 1));

          simulate_imh_proposal_index = simulate_imh_proposal_index + 1;
        } else {
          stop("No imh_proposal found.");
        }
      } else if ((*i) == 3) {
        SymmetricProposalKernel *proposal;

        if (model.containsElementNamed("m_proposal")) {
          List proposal_infos = model["m_proposal"];

          if (Rf_isNewList(proposal_infos[simulate_m_proposal_index])) {
            proposal = get_m_proposal(proposal_infos[simulate_m_proposal_index],
                                      model_parameters, data);
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
          List current_proposal_info = proposal_infos[simulate_m_proposal_index];
          move_factors.push_back(set_move_factor_indices(
              mcmc, current_proposal_info, "m", full_index,
              check_factor_index, simulate_m_proposal_index + 1));

          simulate_m_proposal_index = simulate_m_proposal_index + 1;
        } else {
          stop("No m_proposal found.");
        }
      } else if ((*i) == 4) {
        ProposalKernel *proposal;

        if (model.containsElementNamed("unadjusted_proposal")) {
          List proposal_infos = model["unadjusted_proposal"];
          if (Rf_isNewList(
                  proposal_infos[simulate_unadjusted_proposal_index])) {
            proposal = get_unadjusted_proposal(
                proposal_infos[simulate_unadjusted_proposal_index],
                model_parameters, data);
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
          List current_proposal_info = proposal_infos[simulate_unadjusted_proposal_index];
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

    if (mcmc_weights_method.containsElementNamed("func")) {
      NumericVector mcmc_weights;

      SEXP mcmc_weights_SEXP = mcmc_weights_method["func"];
      mcmc_weights = load_numeric_vector(mcmc_weights_SEXP);

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

std::vector<Parameters> make_initial_points(const List &initial_values) {
  std::vector<Parameters> output;
  output.reserve(initial_values.size());

  for (size_t i = 0; i < initial_values.size(); ++i) {
    Parameters current_parameters;
    if (Rf_isNewList(initial_values[i])) {
      List list_params = initial_values[i];
      if (list_params.size() > 0) {
        CharacterVector names = list_params.names();

        for (size_t i = 0; i < names.size(); ++i) {
          if (isDoubleInList(list_params, Rcpp::as<std::string>(names[i]))) {
            double param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          } else if (isNumericVectorInList(list_params,
                                           Rcpp::as<std::string>(names[i]))) {
            NumericVector param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] =
                as<arma::colvec>(param);
          } else if (isNumericMatrixInList(list_params,
                                           Rcpp::as<std::string>(names[i]))) {
            arma::mat param = list_params[Rcpp::as<std::string>(names[i])];
            current_parameters[Rcpp::as<std::string>(names[i])] = param;
          } else {
            stop("Initial value must be a float/double, vector or matrix.");
          }
        }
      } else {
        Rcpp::stop("List of initial values has empty entries.");
      }

      output.push_back(current_parameters);
    } else {
      Rcpp::stop("List of initial values must be a list (over chains) of "
                 "lists, each of which contain the initial parameters.");
    }
  }

  return output;
}

// [[Rcpp::export]]
size_t ilike_rdtsc() { return rdtsc(); }

Parameters make_algorithm_parameters(const List &algorithm_parameter_list) {

  return list_to_parameters(algorithm_parameter_list);
}

SMCCriterion *get_resampling_method(const List &model,
                                    const List &model_parameters,
                                    size_t number_of_particles) {
  SMCCriterion *smc_method;

  if (model.containsElementNamed("method")) {
    List method_info = model["method"];

    for (size_t i = 0; i < method_info.size(); ++i) {
      List current_method = method_info[i];
      if (current_method.containsElementNamed("adaptive_resampling")) {
        List adaptive_resampling_method = current_method["adaptive_resampling"];

        if (adaptive_resampling_method.containsElementNamed("type") &&
            adaptive_resampling_method.containsElementNamed("parameters")) {
          std::string method = adaptive_resampling_method["type"];
          if (iequals(method, "ess")) {
            if (Rf_isNewList(adaptive_resampling_method["parameters"])) {
              List values = adaptive_resampling_method["parameters"];

              if (values.length() == 1) {
                double proportion =
                    extract_double_parameter(values, model_parameters, 0);
                smc_method = new ESSSMCCriterion(proportion *
                                                 double(number_of_particles));
                return smc_method;
              } else {
                stop("Adaptive resampling using ESS requires specification of "
                     "the proportion of the total number of particles.");
              }
            } else {
              stop("Adaptive resampling using ESS requires specification of "
                   "the proportion of the total number of particles.");
            }
          } else {
            stop("No valid method found for adaptive resampling.");
          }
        }
      }
    }
  }

  inform_once("No method set for adaptive resampling: defaulting to resampling whenever the ESS drops below the number of particles.");
  return new ESSSMCCriterion(double(number_of_particles));
}

SMCCriterion *get_adaptive_target_method(const List &model,
                                         const List &model_parameters,
                                         const List &adaptive_target_method,
                                         size_t number_of_particles) {
  SMCCriterion *smc_method;

  if (adaptive_target_method.containsElementNamed("method")) {
    std::string method = adaptive_target_method["method"];
    if (iequals(method, "ess")) {
      if (adaptive_target_method.containsElementNamed("values")) {
        if (Rf_isNewList(adaptive_target_method["values"])) {
          List values = adaptive_target_method["values"];

          if (values.length() == 1) {
            double proportion =
                extract_double_parameter(values, model_parameters, 0);
            smc_method =
                new ESSSMCCriterion(proportion * double(number_of_particles));
          } else {
            stop("Adaptive target using ESS requires specification of the "
                 "proportion of the total number of particles.");
          }
        } else {
          stop("No valid method found for adaptive target.");
        }
      } else {
        stop("Adaptive target using ESS requires specification of the "
             "proportion of the total number of particles (values need to be "
             "specified).");
      }
    } else if (iequals(method, "cess")) {
      if (adaptive_target_method.containsElementNamed("values")) {
        if (Rf_isNewList(adaptive_target_method["values"])) {
          List values = adaptive_target_method["values"];

          if (values.length() == 1) {
            double proportion =
                extract_double_parameter(values, model_parameters, 0);

            smc_method =
                new CESSSMCCriterion(proportion * double(number_of_particles));
          } else {
            stop("Adaptive target using CESS requires specification of the "
                 "proportion of the total number of particles.");
          }
        } else {
          stop("No valid method found for adaptive target.");
        }
      } else {
        stop("Adaptive target using ESS requires specification the proportion "
             "of the total number of particles (values need to be specified).");
      }
    } else if (iequals(method, "positive")) {
      smc_method = new PositiveSMCCriterion();
    } else {
      stop("No valid method found for adaptive target.");
    }
  } else {
    stop("No valid method found for adaptive target (need method).");
  }

  return smc_method;
}

SMCTermination *get_smc_termination_method(const List &model,
                                           const List &model_parameters,
                                           const List &smc_termination_method) {
  SMCTermination *smc_termination;

  if (smc_termination_method.length() == 0) {
    smc_termination = NULL;
    return smc_termination;
  }

  if (smc_termination_method.containsElementNamed("method")) {
    std::string method = smc_termination_method["method"];
    if (iequals(method, "stable")) {
      if (smc_termination_method.containsElementNamed("values")) {
        if (Rf_isNewList(smc_termination_method["values"])) {
          List values = smc_termination_method["values"];

          if (values.length() == 2) {
            double threshold =
                extract_double_parameter(values, model_parameters, 0);
            int number_in_a_row =
                extract_int_parameter(values, model_parameters, 1);
            if (number_in_a_row < 1)
              stop("ilike::stable: the number of iterations in a row must be at least 1.");
            // The constructor takes (number_in_a_row, threshold).
            smc_termination =
                new StableSMCTermination(size_t(number_in_a_row), threshold);
          } else {
            stop("SMC termination via stability requires the specificaton of a "
                 "threshold and the number of iterations in a row that this "
                 "threshold must be exceeded.");
          }
        } else {
          stop("No valid method found for SMC termination.");
        }
      } else {
        stop("SMC termination via stability requires specification of the "
             "desired ESS (values need to be specified).");
      }
    } else if (iequals(method, "always")) {
      smc_termination = new AlwaysSMCTermination();
    } else {
      stop("No valid method found for SMC termination.");
    }
  } else {
    stop("No valid method found for SMC termination (method needed).");
  }

  return smc_termination;
}

EnsembleShifter *get_enk_shifter_method(const List &model,
                                        const List &model_parameters) {
  EnsembleShifter *enk_shifter;

  if (model.containsElementNamed("method")) {
    List method_info = model["method"];

    for (size_t i = 0; i < method_info.size(); ++i) {
      List current_method = method_info[i];
      if (current_method.containsElementNamed("enk_shifter")) {
        List method = current_method["enk_shifter"];

        if (method.containsElementNamed("type")) {
          std::string type = method["type"];

          if (iequals(type, "stochastic")) {
            enk_shifter = new StochasticEnsembleShifter();
            return enk_shifter;
          } else if (iequals(type, "sqrt")) {
            enk_shifter = new SquareRootEnsembleShifter();
            return enk_shifter;
          } else if (iequals(type, "adjustment")) {
            enk_shifter = new AdjustmentEnsembleShifter();
            return enk_shifter;
          } else {
            stop("No valid method found for EnK shifter.");
          }
        }
      }
    }
  }

  inform_once("No method set for shifting ensemble in EnK: defaulting to stochastic approach.");
  return new StochasticEnsembleShifter();
}

List get_smc_sequencer_info(const List &model, const List &model_parameters,
                            const List &smc_sequencer_method) {
  if (smc_sequencer_method.containsElementNamed("types") &&
      smc_sequencer_method.containsElementNamed("variables") &&
      smc_sequencer_method.containsElementNamed("schedules")) {
    int number_of_bisections;
    if (smc_sequencer_method.containsElementNamed("bisections")) {
      List values = smc_sequencer_method["bisections"];
      number_of_bisections = extract_int_parameter(values, model_parameters, 0);
    } else {
      inform_once("Bisections not set in adaptive target method, using 100.");
      number_of_bisections = 100;
    }

    StringVector types;
    if (isStringVectorInList(smc_sequencer_method, "types")) {
      types = smc_sequencer_method["types"];
    }

    StringVector variables;
    if (isStringVectorInList(smc_sequencer_method, "variables")) {
      variables = smc_sequencer_method["variables"];
    }

    List schedules;
    if (Rf_isNewList(smc_sequencer_method["schedules"])) {
      schedules = smc_sequencer_method["schedules"];
    } else if (isNumericVectorInList(smc_sequencer_method, "schedules")) {
      schedules = List::create(smc_sequencer_method["schedules"]);
    } else {
      stop("Invalid specification of schedules for SMC: must be either a "
           "numeric vector or a list of numeric vectors.");
    }

    if ((types.size() == variables.size()) &&
        (types.size() == schedules.size())) {
      std::vector<std::string> types_for_output;
      types_for_output.reserve(types.size());
      for (size_t i = 0; i < types.size(); ++i) {
        types_for_output.push_back(String(types[i]).get_cstring());
      }

      std::vector<std::string> variables_for_output;
      variables_for_output.reserve(variables.size());
      for (size_t i = 0; i < variables.size(); ++i) {
        variables_for_output.push_back(String(variables[i]).get_cstring());
      }

      std::vector<std::vector<double>> schedules_for_output;
      schedules_for_output.reserve(schedules.size());
      for (size_t i = 0; i < schedules.size(); ++i) {
        std::vector<double> schedule_for_output;

        if (isNumericVectorInList(schedules, i + 1)) {
          NumericVector schedule = schedules[i];
          schedule_for_output = as<std::vector<double>>(schedule);
        } else {
          stop("Schedule invalid: needs to be a numeric vector.");
        }

        schedules_for_output.push_back(schedule_for_output);
      }

      if (smc_sequencer_method.containsElementNamed("factors_affected")) {
        if (isNumericVectorInList(smc_sequencer_method, "factors_affected")) {
          NumericVector factors_affected =
              smc_sequencer_method["factors_affected"];

          return List::create(types_for_output, variables_for_output,
                              schedules_for_output, number_of_bisections,
                              factors_affected);
        } else {
          stop("Factors affected invalid: needs to be a numeric vector.");
        }
      } else {
        return List::create(types_for_output, variables_for_output,
                            schedules_for_output, number_of_bisections);
      }
    } else {
      stop("Invalid specification of schedules for SMC: types, variables and "
           "schedules (list) must all be the same length.");
    }
  } else {
    std::vector<std::string> types_for_output;
    std::vector<std::string> variables_for_output;
    std::vector<std::vector<double>> schedules_for_output;
    size_t number_of_bisections = 100;
    return List::create(types_for_output, variables_for_output,
                        schedules_for_output, number_of_bisections);
    // stop("No valid method found for SMC sequence: need to specify types,
    // variables and schedules.");
  }
}
