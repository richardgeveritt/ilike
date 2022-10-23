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
//#include "ensemble_member.h"

class ModelAndAlgorithm;
class LikelihoodEstimatorOutput;
class Transform;
class ProposalKernel;
class GradientEstimatorOutput;
class GradientEstimator;
class FactorVariables;
class EnsembleFactorVariables;

class Particle
{

public:

  Particle();
  Particle(const Parameters &parameters_in);
  Particle(const Parameters &parameters_in,
           FactorVariables* factor_variables_in);
  Particle(const Parameters &parameters_in,
           FactorVariables* factor_variables_in,
           double previous_target_evaluated_in);
  Particle(const Parameters &parameters_in,
           EnsembleFactorVariables* ensemble_factor_variables_in);
  virtual ~Particle();

  Particle(const Particle &another);
  void operator=(const Particle &another);

  friend std::ostream& operator<<(std::ostream& os, const Particle &p);
  
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
  void evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                             const Parameters &conditioned_on_parameters);
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                              const Parameters &conditioned_on_parameters);
  double evaluate_likelihoods(const Index* index);
  double evaluate_likelihoods(const Index* index,
                              const Parameters &conditioned_on_parameters);
  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index,
                                                       const Parameters &conditioned_on_parameters);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index,
                                                                        const Parameters &conditioned_on_parameters);
  double subsample_evaluate_likelihoods(const Index* index);
  double subsample_evaluate_likelihoods(const Index* index,
                                        const Parameters &conditioned_on_parameters);
  
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double incremental_temperature);
  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double incremental_temperature,
                                             const Parameters &conditioned_on_parameters);
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double incremental_temperature,
                                                       const Parameters &conditioned_on_parameters);
  
  double evaluate_all_likelihoods(const Index* index);
  double evaluate_all_likelihoods(const Index* index,
                                  const Parameters &conditioned_on_parameters);
  double subsample_evaluate_all_likelihoods(const Index* index,
                                            const Parameters &conditioned_on_parameters);
  
  arma::colvec get_vector() const;
  arma::colvec get_vector(const std::vector<std::string> &state_names) const;
  
  void set_move_transformed_parameters();
  void set_move_transformed_parameters(Transform* transform_in);
  
  Parameters parameters;
  
  // parameters after they've been put through the transformation of the move recent move
  Parameters move_transformed_parameters;
  
  // not stored here - this points to one of the above, depending on which one if being used by the current move
  Parameters* move_parameters;
  
  // Parameters stored here, proposal kernel is pointed to.
  boost::unordered_map< const ProposalKernel*, GradientEstimatorOutput*> gradient_estimator_outputs;
  
  // pointer returned will be owned by this class
  // function is member of Variables since we will not always create a new GradientEstimatorOutput if one already exists for this proposal
  GradientEstimatorOutput* initialise_gradient_estimator_output(const ProposalKernel* proposal,
                                                                GradientEstimator* gradient_estimator);
  
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
                                       const Index* index);
  arma::mat direct_get_gradient_of_log(const std::string &variable,
                                       const Index* index,
                                       const Parameters &conditioned_on_parameters);
  
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                 const Index* index);
  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                 const Index* index,
                                                 const Parameters &conditioned_on_parameters);
  
  void simulate_factor_variables(Particle* previous_particle);
  
  void simulate_ensemble_factor_variables(Particle* previous_particle);
  
  void simulate_factor_variables(Particle* previous_particle,
                                 const Parameters &conditioned_on_parameters);
  
  void simulate_ensemble_factor_variables(Particle* previous_particle,
                                          const Parameters &conditioned_on_parameters);
  
  // not stored here
  Particle* previous_self;

protected:
  
  // not stored here
  Transform* move_transform;

  void make_copy(const Particle &another);

};

#endif
