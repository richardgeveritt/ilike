#ifndef PARTICLES_H
#define PARTICLES_H

//#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
//using namespace RcppParallel;

#include "particle.h"
#include "distributions.h"

namespace ilike
{
class SMC;
class SequentialSMCWorker;
class MoveOutput;
class Factors;

class Particles
{
public:
  
  Particles();
  Particles(size_t number_of_particles_in);
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
  
  virtual ~Particles();
  
  Particles(const Particles &another);
  Particles& operator=(const Particles &another);
  
  Particles(Particles &&another);
  Particles& operator=(Particles &&another);
  
  void reserve(size_t number_of_particles_in);
  //void push_back(const Parameters &parameters_in,
  //               Factors* factors_in);
  //void push_back(const Particle &particle_in);
  void push_back(MoveOutput* move_output_in);
  void push_back(Parameters &&parameters_in,
                 Factors* factors_in,
                 const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                 const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  void push_back(Particle &&particle_in);
  void push_back(const Parameters &parameters_in,
                 Factors* factors_in,
                 const std::vector<const ProposalKernel*>* proposals_to_transform_for_in,
                 const std::vector<const ProposalKernel*>* proposals_to_find_gradient_for_in);
  void push_back(const Particle &particle_in);
  
  Particle* add_particle();
  
  void simulate_resampling_variables(RandomNumberGenerator &rng);
  
  //arma::colvec get_resampling_variables() const;
  //void set_ancestor_variables(const std::vector<size_t> &ancestor_variables_in);
  
  void initialise_weights();
  void update_weights(const arma::colvec &latest_unnormalised_log_incremental_weights);
  double calculate_log_normalising_constant();
  void normalise_weights();
  void resample();
  
  arma::rowvec get_output_lengths() const;
  
  void set_previous_target_evaluated_to_target_evaluated();
  void subsample_set_previous_target_evaluated_to_target_evaluated();
  
  arma::mat get_most_recent_matrix_particles(const std::vector<std::string> &vector_variables) const;
  
  //virtual List& operator[](int index)=0;
  //virtual List operator[](int index) const=0;
  
  friend std::ostream& operator<<(std::ostream& os, const Particles &p);
  
  // Access only elements in std::vector<Particle>.
  //double& operator[](const size_t &i);
  //double operator[](const size_t &i) const;
  
  MoveOutput* operator[](const size_t &i);
  MoveOutput* operator[](const size_t &i) const;
  
  MoveOutput* back();
  MoveOutput* back() const;
  
  size_t size() const;
  
  // Outer vector is over particles, inner is stored here
  std::vector< MoveOutput* > particles;
  
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
  
  double ess;
  
  void close_ofstreams();
  
protected:
  
  //friend SequentialSMCWorker;
  
  void make_copy(const Particles &another);
  
  void make_copy(Particles &&another);
  
  // Stored here. // moved to model_and_algorithm
  //std::vector<LikelihoodEstimator*> likelihood_estimators;
  
};
}

#endif
