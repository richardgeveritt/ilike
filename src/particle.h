#ifndef PARTICLE_H
#define PARTICLE_H

//#include <Rcpp.h>

#include <RcppArmadillo.h>
using namespace Rcpp;
#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>
#include "parameters.h"
//#include "variables.h"
#include "index.h"
#include "distributions.h"
#include "proposal_store.h"

class ModelAndAlgorithm;
class LikelihoodEstimatorOutput;
class Transform;
class ProposalKernel;
class GradientEstimatorOutput;
class GradientEstimator;
class Factors;
class EnsembleFactors;
class FactorVariables;
class EnsembleFactorVariables;

class Particle
{

public:

  Particle();
  //Particle(const Parameters &&parameters_in);
  
  Particle(const Parameters &parameters_in);
  
  Particle(const Parameters &parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  
  Particle(const Parameters &parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
           const Parameters &conditioned_on_parameters);
  
  Particle(const Parameters &parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);
  
  Particle(const Parameters &parameters_in,
           FactorVariables* factor_variables_in);
  
  Particle(const Parameters &parameters_in,
           FactorVariables* factor_variables_in,
           double previous_target_evaluated_in);
  
  Particle(const Parameters &parameters_in,
           const EnsembleFactors* ensemble_factors_in);
  
  Particle(const Parameters &parameters_in,
           const EnsembleFactors* ensemble_factors_in,
           const Parameters &conditioned_on_parameters);
  
  Particle(const Parameters &parameters_in,
           const EnsembleFactors* ensemble_factors_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);
  
  Particle(const Parameters &parameters_in,
           EnsembleFactorVariables* ensemble_factor_variables_in);
  
  
  Particle(Parameters &&parameters_in);
  
  Particle(Parameters &&parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  
  Particle(Parameters &&parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
           const Parameters &conditioned_on_parameters);
  
  Particle(Parameters &&parameters_in,
           const Factors* factors_in,
           const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
           const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);
  
  Particle(Parameters &&parameters_in,
           FactorVariables* factor_variables_in);
  
  Particle(Parameters &&parameters_in,
           FactorVariables* factor_variables_in,
           double previous_target_evaluated_in);
  
  Particle(Parameters &&parameters_in,
           const EnsembleFactors* ensemble_factors_in);
  
  Particle(Parameters &&parameters_in,
           const EnsembleFactors* ensemble_factors_in,
           const Parameters &conditioned_on_parameters);
  
  Particle(Parameters &&parameters_in,
           const EnsembleFactors* ensemble_factors_in,
           const Parameters &conditioned_on_parameters,
           const Parameters &sequencer_parameters);
  
  Particle(Parameters &&parameters_in,
           EnsembleFactorVariables* ensemble_factor_variables_in);
  
  
  virtual ~Particle();
  
  void setup(Factors* factors_in);
  void setup(Factors* factors_in,
             const Parameters &conditioned_on_parameters);
  void setup(Factors* factors_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  void setup(const Parameters &parameters_in,
             Factors* factors_in);
  void setup(const Parameters &parameters_in,
             Factors* factors_in,
             const Parameters &conditioned_on_parameters);
  void setup(const Parameters &parameters_in,
             Factors* factors_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  void setup(const Parameters &parameters_in,
             EnsembleFactors* factors_in);
  void setup(const Parameters &parameters_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters);
  void setup(const Parameters &parameters_in,
             EnsembleFactors* factors_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  
  Particle copy_without_factor_variables() const;

  Particle(const Particle &another);
  Particle& operator=(const Particle &another);
  
  Particle(Particle &&another);
  Particle& operator=(Particle &&another);

  friend std::ostream& operator<<(std::ostream& os, const Particle &p);
  
  void simulate_factor_variables();
  void simulate_ensemble_factor_variables();
  void simulate_proposal_variables(const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  
  /*
  void evaluate_smcfixed_part_of_likelihoods();
  void evaluate_smcfixed_part_of_likelihoods(const Parameters &conditioned_on_parameters);
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods();
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Parameters &conditioned_on_parameters);
  double evaluate_likelihoods();
  double evaluate_likelihoods(const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Parameters &conditioned_on_parameters);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Parameters &conditioned_on_parameters);
  double subsample_evaluate_likelihoods();
  double subsample_evaluate_likelihoods(const Parameters &conditioned_on_parameters);
  */
  
  void evaluate_smcfixed_part_of_likelihoods(const Index* index);
  /*
  void evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                             const Parameters &conditioned_on_parameters);
  */
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  
  /*
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                              const Parameters &conditioned_on_parameters);
  */
  
  double evaluate_likelihoods(const Index* index);
  
  /*
  double evaluate_likelihoods(const Index* index,
                              const Parameters &conditioned_on_parameters);
  */
  
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  
  /*
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                       const Parameters &conditioned_on_parameters);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                        const Parameters &conditioned_on_parameters);
  */
  
  double subsample_evaluate_likelihoods(const Index* index);
  
  /*
  double subsample_evaluate_likelihoods(const Index* index,
                                        const Parameters &conditioned_on_parameters);
  */
  
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double incremental_temperature);
  /*
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double incremental_temperature,
                                             const Parameters &conditioned_on_parameters);
  */
  
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double incremental_temperature);
  
  /*
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double incremental_temperature,
                                                       const Parameters &conditioned_on_parameters);
  */
  
  double evaluate_all_likelihoods(const Index* index);
  
  /*
  double evaluate_all_likelihoods(const Index* index,
                                  const Parameters &conditioned_on_parameters);
  */
  
  double subsample_evaluate_all_likelihoods(const Index* index);
  
  /*
  double subsample_evaluate_all_likelihoods(const Index* index,
                                            const Parameters &conditioned_on_parameters);
  */
  
  //arma::colvec get_vector() const;
  
  arma::colvec get_colvec(const std::string &state_name) const;
  
  arma::rowvec get_rowvec(const std::string &state_name) const;
  
  arma::colvec get_colvec(const std::vector<std::string> &state_names) const;
  
  arma::rowvec get_rowvec(const std::vector<std::string> &state_names) const;
  
  Parameters get_transformed_parameters(const ProposalKernel* proposal_in) const;
  
  GradientEstimatorOutput* get_gradient_estimator_output(const ProposalKernel* proposal_in) const;
  /*
  void set_current_transformed_parameters(const ProposalKernel* proposal_in);
  void set_previous_transformed_parameters(const ProposalKernel* proposal_in,
                                                Transform* transform_in);
  
  void set_previous(const Particle &previous_particle);
  
  Parameters get_previous_transformed_parameters(const ProposalKernel* proposal_in) const;
  */
   
  //void set_previous_move_transformed_parameters();
  //void set_previous_move_transformed_parameters(Transform* transform_in);
  
  Parameters parameters;
  
  // parameters after they've been put through the transformation of the move recent move
  //Parameters move_transformed_parameters;
  //Parameters previous_move_transformed_parameters;
  
  boost::unordered_map< const ProposalKernel*, ProposalStore> current_proposal_store;
  
  const std::vector<const ProposalKernel*>* proposals_to_transform_for_pointer;
  const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_pointer;
  
  //boost::unordered_map< const ProposalKernel*, ProposalStore> previous_proposal_store;
  
  // not stored here - this points to one of the above, depending on which one if being used by the current move
  //Parameters* move_parameters;
  
  // not stored here
  //Transform* move_transform;
  
  // Parameters stored here, proposal kernel is pointed to.
  //boost::unordered_map< const ProposalKernel*, GradientEstimatorOutput*> gradient_estimator_outputs;
  
  // pointer returned will be owned by this class
  // function is member of Variables since we will not always create a new GradientEstimatorOutput if one already exists for this proposal
  /*
  GradientEstimatorOutput* initialise_gradient_estimator_output(const ProposalKernel* proposal,
                                                                GradientEstimator* gradient_estimator);
  
  GradientEstimatorOutput* initialise_previous_gradient_estimator_output(const ProposalKernel* proposal,
                                                                         GradientEstimator* gradient_estimator);
  */
  
  // take the transformed parameters, store them in the ProposalInfo, and also go through the inverse transform to set the parameters
  //void set_parameters(const ProposalKernel* proposal,
  //                    const Parameters &proposed_transformed);
  
  // Proposal kernel is pointed to
  boost::unordered_map< const ProposalKernel*, bool> accepted_outputs;
  
  void set_acceptance(const ProposalKernel* proposal_in,
                      bool accepted_in);
  
  void erase_mcmc_adaptation_info();
  
  // What happened at the last target.
  double previous_target_evaluated;
  double subsample_previous_target_evaluated;
  
  // What happened at the current step (whether the true target or the one we are trying in adaptive SMC).
  double target_evaluated;
  double subsample_target_evaluated;
  
  double previous_ensemble_target_evaluated;
  double subsample_previous_ensemble_target_evaluated;
  
  double ensemble_target_evaluated;
  double subsample_ensemble_target_evaluated;
  
  // Stored here.
  FactorVariables* factor_variables;
  EnsembleFactorVariables* ensemble_factor_variables;
  
  arma::mat direct_get_gradient_of_log(const std::string &variable,
                                       const Index* index) const;
  
  /*
  arma::mat direct_get_gradient_of_log(const std::string &variable,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters);
  */
  
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                 const Index* index) const;
  
  /*
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                 const Index* index,
                                                 const Parameters &conditioned_on_parameters);
  */
  
  void simulate_factor_variables(const Factors* factors);
  
  void simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors);
  
  /*
  void simulate_factor_variables(Factors* factors,
                                 const Parameters &conditioned_on_parameters);
  
  void simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors,
                                          const Parameters &conditioned_on_parameters);
  */
  
  void subsample_simulate_factor_variables(const Factors* factors);
  
  void subsample_simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors);
  
  /*
  void subsample_simulate_factor_variables(Factors* factors,
                                           const Parameters &conditioned_on_parameters);
  
  void subsample_simulate_ensemble_factor_variables(EnsembleFactors* ensemble_factors,
                                                    const Parameters &conditioned_on_parameters);
  */
  
  // not stored here
  const Particle* previous_self;

protected:

  void make_copy(const Particle &another);
  
  void make_copy(Particle &&another);

};

#endif
