#ifndef PARTICLES_H
#define PARTICLES_H

//#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
//using namespace RcppParallel;

#include "particle.h"
#include "sample_store.h"
#include "distributions.h"

namespace ilike
{
  /**
   * @file particles.h
   * @brief Defines the SMC class.
   *
   * Provides smc functionality.
   *
   * @namespace ilike
   * @class SMC
   * @brief The smc class.
   */


class SMC;
class SequentialSMCWorker;
class MoveOutput;
class Factors;
class FactorVariables;
class Index;

class Particles
{
public:
  
  /**
   * @brief Performs the particles operation.
   */
  Particles();
  /**
   * @brief Performs the particles operation.
   *
   * @param number_of_particles_in The number of particles.
   */
  Particles(size_t number_of_particles_in);
  /**
   * @brief Performs the particles operation.
   *
   * @param particles_in The particles.
   */
  Particles(const std::vector< MoveOutput* > &particles_in);
  /*
   Particles(std::vector<Parameters> &initial_values_in,
   const arma::colvec &log_probabilities_of_initial_values_in,
   Factors* factors_in,
   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
   */
  Particles(std::vector<Parameters> &initial_values_in,
            const arma::colvec &log_probabilities_of_initial_values_in,
            Factors* factors_in,
            const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
            const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
            const Parameters &conditioned_on_parameters);
  Particles(std::vector<Parameters> &initial_values_in,
            const arma::colvec &log_probabilities_of_initial_values_in,
            Factors* factors_in,
            const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
            const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
            const Parameters &conditioned_on_parameters,
            const Parameters &sequencer_parameters);
  
  /*
   void setup(std::vector<Parameters> &initial_values_in,
   const arma::colvec &log_probabilities_of_initial_values_in,
   Factors* factors_in,
   const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
   const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
   */
  
  void setup(std::vector<Parameters> &initial_values_in,
             const arma::colvec &log_probabilities_of_initial_values_in,
             Factors* factors_in,
             const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
             const Parameters &conditioned_on_parameters);
  
  void setup(std::vector<Parameters> &initial_values_in,
             const arma::colvec &log_probabilities_of_initial_values_in,
             Factors* factors_in,
             const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
             const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in,
             const Parameters &conditioned_on_parameters,
             const Parameters &sequencer_parameters);
  
  /**
   * @brief Performs the ~particles operation.
   */
  virtual ~Particles();
  
  /**
   * @brief Performs the particles operation.
   *
   * @param another The SMC instance to copy from.
   */
  Particles(const Particles &another);
  /**
   * @brief Assignment operator for SMC.
   *
   * @param another The SMC instance to copy from.
   */
  Particles& operator=(const Particles &another);
  
  /**
   * @brief Performs the particles operation.
   *
   * @param another The SMC instance to copy from.
   */
  Particles(Particles &&another);
  /**
   * @brief Assignment operator for SMC.
   *
   * @param another The SMC instance to copy from.
   */
  Particles& operator=(Particles &&another);
  
  /**
   * @brief Performs the reserve operation.
   *
   * @param number_of_particles_in The number of particles.
   */
  void reserve(size_t number_of_particles_in);
  //void push_back(const Parameters &parameters_in,
  //               Factors* factors_in);
  //void push_back(const Particle &particle_in);
  /**
   * @brief Performs the push back operation.
   *
   * @param move_output_in The move output.
   */
  void push_back(MoveOutput* move_output_in);
  void push_back(Parameters &&parameters_in,
                 Factors* factors_in,
                 const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                 const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  /**
   * @brief Performs the push back operation.
   *
   * @param particle_in The particle.
   */
  void push_back(Particle &&particle_in);
  void push_back(const Parameters &parameters_in,
                 Factors* factors_in,
                 const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                 const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  /**
   * @brief Performs the push back operation.
   *
   * @param particle_in The particle.
   */
  void push_back(const Particle &particle_in);
  
  /**
   * @brief Performs the add particle operation.
   *
   * @return The result.
   */
  Particle* add_particle();
  
  /**
   * @brief Simulates resampling variables.
   *
   * @param rng The rng.
   */
  void simulate_resampling_variables(RandomNumberGenerator &rng);
  
  //arma::colvec get_resampling_variables() const;
  //void set_ancestor_variables(const std::vector<size_t> &ancestor_variables_in);
  
  /**
   * @brief Performs the initialise weights operation.
   */
  void initialise_weights();
  /**
   * @brief Updates weights.
   *
   * @param latest_unnormalised_log_incremental_weights The latest unnormalised log incremental weights.
   */
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  /**
   * @brief Performs the calculate log normalising constant operation.
   *
   * @return The result.
   */
  double calculate_log_normalising_constant();
  /**
   * @brief Performs the normalise weights operation.
   */
  void normalise_weights();
  /**
   * @brief Resamples the particle population.
   */
  void resample();
  
  /**
   * @brief Returns the output lengths.
   *
   * @return The result.
   */
  arma::rowvec get_output_lengths() const;
  
  /**
   * @brief Sets the previous target evaluated to target evaluated.
   */
  void set_previous_target_evaluated_to_target_evaluated();
  /**
   * @brief Performs the subsample set previous target evaluated to target evaluated operation.
   */
  void subsample_set_previous_target_evaluated_to_target_evaluated();
  
  /**
   * @brief Returns the most recent matrix particles.
   *
   * @param vector_variables The vector variables.
   *
   * @return The result.
   */
  arma::mat get_most_recent_matrix_particles(const std::vector<std::string> &vector_variables) const;
  
  //virtual List& operator[](int index)=0;
  //virtual List operator[](int index) const=0;
  
  /**
   * @brief Performs the operator<< operation.
   *
   * @param os The os.
   * @param p The p.
   *
   * @return The result.
   */
  friend std::ostream& operator<<(std::ostream& os, const Particles &p);
  
  // Access only elements in std::vector<Particle>.
  //double& operator[](const size_t &i);
  //double operator[](const size_t &i) const;
  
  /**
   * @brief Performs the operator[] operation.
   *
   * @param i The i.
   *
   * @return The result.
   */
  MoveOutput* operator[](const size_t &i);
  /**
   * @brief Performs the operator[] operation.
   *
   * @param i The i.
   *
   * @return The result.
   */
  MoveOutput* operator[](const size_t &i) const;
  
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  MoveOutput* back();
  /**
   * @brief Performs the back operation.
   *
   * @return The result.
   */
  MoveOutput* back() const;
  
  /**
   * @brief Performs the size operation.
   *
   * @return The result.
   */
  size_t size() const;
  
  // Outer vector is over particles, inner is stored here
  std::vector< MoveOutput* > particles;

  // Objects of a dropped generation, given to this one to be reused for its
  // new particles (dev/plan.md, Stage 3 (a); used by the particle filter).
  // Owned: deleted with this object if not taken. Not copied.
  std::vector< MoveOutput* > spare;
  MoveOutput* take_spare();
  
  arma::colvec unnormalised_log_weights;
  
  arma::colvec incremental_log_weights;
  
  arma::colvec normalised_log_weights;
  
  arma::colvec previous_normalised_log_weights;
  
  double log_normalising_constant_ratio;
  
  // Also store anything that needs to be stored in a different structure across all particles.
  arma::colvec resampling_variables;
  std::vector<size_t> ancestor_variables;
  
  bool resampled_flag;
  
  // only here for output
  Parameters schedule_parameters;
  
  // NaN until computed (MCMC output used to write uninitialised memory here).
  double ess = arma::datum::nan;

  // Flat form of a generation (dev/plan.md 2.6a): the points kept as numbers
  // instead of one MoveOutput (a rich Particle) each. Used by algorithms that
  // simulate and weight a generation in chunks with reused particle objects
  // (chunked importance sampling). In the flat form `particles` is empty, and
  // the accessors that hand out a rich particle (operator[], back()) stop, so
  // that a reader the flat form does not serve is found by a message.
  //  - flat_points: the variables of every point (and, as extra columns, the
  //    compact factor records, when the factor output is to be written and
  //    every factor has one);
  //  - flat_factor_variables: when the factor output is to be written but a
  //    factor has no compact record, each point's factor variables, kept whole
  //    (owned; detached from any particle);
  //  - flat_record_writer: a factor-variables object of the same factors, used
  //    only to write the compact records (owned);
  //  - flat_incremental_log_weights: each point's incremental log weight at
  //    the generation's target (what weight_for_adapting_sequence returns).
  bool flat = false;
  size_t flat_size = 0;
  SampleStore flat_points;
  std::vector<FactorVariables*> flat_factor_variables;
  FactorVariables* flat_record_writer = NULL;
  arma::colvec flat_incremental_log_weights;
  //  - for a generation made by MCMC moves (dev/plan.md 2.6a, G3):
  //    flat_points holds every stored chain state, particle by particle,
  //    flat_output_lengths the number stored for each particle, and
  //    flat_final_points the final state of each particle's chain (from
  //    which the next moves start, with the adaptive summaries below).
  std::vector<double> flat_output_lengths;
  //  - flat_records_from_summaries: the compact factor records were taken
  //    before the generation's target was chosen (a first SMC generation),
  //    and are recomputed from the summaries once it has been
  //    (SMC::finish_weighting_from_summaries), as the per-particle path wrote
  //    the values its last weighting left.
  bool flat_records_from_summaries = false;
  SampleStore flat_final_points;
  bool is_flat() const { return this->flat; }
  void clear_flat();

  // Adaptive summaries of this generation's points (dev/plan.md 2.6a, G2):
  // each point's adaptive summary over an index (one column per point) and
  // its previous_target_evaluated, so that weights at a trial target are a
  // loop over these numbers (SMC::weight_from_summaries_if_possible).
  // Extracted once, on first use, from the rich particles (whose factor
  // values must be up to date); for a flat generation, filled while it is
  // simulated. Derived data: not copied (a copy extracts them again).
  //  - adaptive_summary_state: 0 not tried, 1 available, -1 not available;
  //  - adaptive_summary_evaluator: factor variables of the same factors used
  //    to evaluate summaries (one of the particles', or flat_summary_evaluator);
  //  - adaptive_summaries_used / last_trial_schedule: weights were last
  //    computed from summaries, at these schedule parameters (so that the
  //    particles' cached values can be brought to that state afterwards).
  arma::mat adaptive_summaries;
  arma::colvec adaptive_previous_targets;
  const Index* adaptive_summary_index = NULL;
  int adaptive_summary_state = 0;
  const FactorVariables* adaptive_summary_evaluator = NULL;
  FactorVariables* flat_summary_evaluator = NULL; // owned
  bool adaptive_summaries_used = false;
  Parameters last_trial_schedule;
  bool prepare_adaptive_summaries(const Index* index);

  // Make the factor-variables objects of a flat generation belong to
  // `factors` (G5); false if the generation holds rich particles.
  bool repoint_flat_to(const Factors* factors,
                       const Index* index = NULL);
  void forget_adaptive_summaries();
  
  /**
   * @brief Closes any open file streams.
   */
  void close_ofstreams();
  
protected:
  
  //friend SequentialSMCWorker;
  
  /**
   * @brief Copies the state of another SMC into this object.
   *
   * @param another The SMC instance to copy from.
   */
  void make_copy(const Particles &another);
  
  /**
   * @brief Copies the state of another SMC into this object.
   *
   * @param another The SMC instance to copy from.
   */
  void make_copy(Particles &&another);
  
  // Stored here. // moved to model_and_algorithm
  //std::vector<LikelihoodEstimator*> likelihood_estimators;
  
};
}

#endif
