#ifndef PARTICLE_H
#define PARTICLE_H

#include <RcppArmadillo.h>
#include <memory>
using namespace Rcpp;
#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>
#include "parameters.h"
#include "index.h"
#include "distributions.h"
#include "proposal_store.h"

namespace ilike
{

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

// Whether the last proposal of each proposal kernel was accepted, keyed by the
// kernel's instance index (read by MCMC adaptors). A particle meets a handful
// of kernels at most, so a vector searched linearly replaces the hash map
// used before: copying an empty record allocates nothing, and a particle is
// copied at every proposal (dev/plan.md, Stage 3 item 2).
class AcceptanceRecords
{
public:
  typedef std::vector<std::pair<int, bool>>::const_iterator const_iterator;

  bool& operator[](int key)
  {
    for (auto i = this->records.begin(); i != this->records.end(); ++i)
      if (i->first == key)
        return i->second;
    this->records.emplace_back(key, false);
    return this->records.back().second;
  }

  const_iterator find(int key) const
  {
    for (auto i = this->records.begin(); i != this->records.end(); ++i)
      if (i->first == key)
        return i;
    return this->records.end();
  }

  const_iterator begin() const { return this->records.begin(); }
  const_iterator end() const { return this->records.end(); }
  size_t size() const { return this->records.size(); }
  void clear() { this->records.clear(); }

private:
  std::vector<std::pair<int, bool>> records;
};

class Particle
{
  
public:
  
  Particle();
  
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
  // Make this particle what Particle(std::move(parameters_in), factors_in,
  // proposals_to_transform_for_in, proposals_to_find_gradient_for_in,
  // conditioned_on_parameters) would be, reusing its factor-variable and
  // output objects where they support it (dev/plan.md 2.6a: the reused
  // particle objects of a chunk). The outputs receive exactly the calls a new
  // object's would.
  void refill(Parameters &&parameters_in,
              const Factors* factors_in,
              const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
              const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
              const Parameters &conditioned_on_parameters,
              const Parameters* second_conditioned_on_parameters = NULL);

  void simulate_proposal_variables(const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                                   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);

  void evaluate_smcfixed_part_of_likelihoods(const Index* index);
  double evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);

  double evaluate_likelihoods(const Index* index);

  // Sum of the cached factor values over index (no evaluation; see
  // FactorVariables::get_cached_likelihoods).
  double get_cached_likelihoods(const Index* index) const;

  void subsample_evaluate_smcfixed_part_of_likelihoods(const Index* index);
  double subsample_evaluate_smcadaptive_part_given_smcfixed_likelihoods(const Index* index);

  double subsample_evaluate_likelihoods(const Index* index);

  double evaluate_ensemble_likelihood_ratios(const Index* index,
                                             double incremental_temperature,
                                             const std::vector<arma::mat> &inv_sigma_precomps,
                                             const std::vector<double> &log_det_precomps) const;
  
  double subsample_evaluate_ensemble_likelihood_ratios(const Index* index,
                                                       double incremental_temperature,
                                                       const std::vector<arma::mat> &inv_sigma_precomps,
                                                       const std::vector<double> &log_det_precomps) const;

  double evaluate_all_likelihoods(const Index* index);

  double subsample_evaluate_all_likelihoods(const Index* index);

  arma::colvec get_colvec(const std::string &state_name) const;
  
  arma::rowvec get_rowvec(const std::string &state_name) const;
  
  arma::colvec get_colvec(const std::vector<std::string> &state_names) const;
  
  arma::rowvec get_rowvec(const std::vector<std::string> &state_names) const;
  
  // The particle's parameters on the scale the proposal works on: the
  // parameters themselves for a proposal without a transform (by reference:
  // returning a copy made every density evaluation and simulation copy the
  // parameter hash map), otherwise the transformed parameters stored for the
  // proposal (computed and kept on first use).
  const Parameters& get_transformed_parameters(const ProposalKernel* proposal_in) const;

  // Exchange the whole state of two particles (every member), keeping each
  // one's factor variables pointing back to the particle that now holds them.
  // Used to accept a proposal held in a scratch particle in O(1).
  void swap(Particle &another);
  
  GradientEstimatorOutput* get_gradient_estimator_output(const ProposalKernel* proposal_in) const;

  Parameters parameters;
  
  // parameters after they've been put through the transformation of the move recent move
  //Parameters move_transformed_parameters;
  //Parameters previous_move_transformed_parameters;
  
  // mutable: get_transformed_parameters (const) keeps a transform it computes
  
  mutable boost::unordered_map< int, ProposalStore> current_proposal_store;
  
  // NULL unless set: several constructors left these uninitialised, and code
  // tests them against NULL before following them (B26)
  const std::vector<const ProposalKernel*>* proposals_to_transform_for_pointer = NULL;
  const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_pointer = NULL;

  // not stored here - this points to one of the above, depending on which one if being used by the current move
  //Parameters* move_parameters;
  
  // not stored here
  //Transform* move_transform;
  
  // Parameters stored here, proposal kernel is pointed to.
  //boost::unordered_map< const ProposalKernel*, GradientEstimatorOutput*> gradient_estimator_outputs;
  
  // pointer returned will be owned by this class
  // function is member of Variables since we will not always create a new GradientEstimatorOutput if one already exists for this proposal
  
  // take the transformed parameters, store them in the ProposalInfo, and also go through the inverse transform to set the parameters
  //void set_parameters(const ProposalKernel* proposal,
  //                    const Parameters &proposed_transformed);
  
  // Proposal kernel is pointed to
  AcceptanceRecords accepted_outputs;
  
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
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<FactorVariables> factor_variables;
  // owned (dev/plan.md 3.5, CH3)
  std::unique_ptr<EnsembleFactorVariables> ensemble_factor_variables;
  
  arma::mat direct_get_gradient_of_log(const std::string &variable,
                                       const Index* index) const;

  arma::mat direct_subsample_get_gradient_of_log(const std::string &variable,
                                                 const Index* index) const;

  void simulate_factor_variables(const Factors* factors);

  // Replace this particle's factor variables with
  // current->copy_with_new_outputs(this->parameters, new_outputs_index).
  void simulate_factor_variables(const FactorVariables* current,
                                 const Index* new_outputs_index);
  
  void simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors);

  void subsample_simulate_factor_variables(const Factors* factors);
  
  void subsample_simulate_ensemble_factor_variables(const EnsembleFactors* ensemble_factors);

  // not stored here
  const Particle* previous_self = NULL;
  
  void tell_factors_to_forget_they_were_already_written_to_file();
  
protected:
  
  void make_copy(const Particle &another);
  
  void make_copy(Particle &&another);
  
};
}

#endif
